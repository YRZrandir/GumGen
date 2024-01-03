#include <chrono>
#include <iostream>
#include <memory>
#include <thread>
//#include <emscripten.h>
#include "gumgen/FakeGumGen.h"
#include <nlohmann/json.hpp>
#include <CGAL/boost/graph/io.h>
#ifdef __EMSCRIPTEN__
#include <pthread.h>
#endif
#ifndef EMSCRIPTEN_KEEPALIVE
#define EMSCRIPTEN_KEEPALIVE
#endif
extern "C"{
std::unique_ptr<FakeGumGen> gumgen = nullptr;
std::unique_ptr<FakeGumGen> gumgen2 = nullptr;

int EMSCRIPTEN_KEEPALIVE CreateGumInit(
    const float* vertices, const unsigned int* indices, unsigned nb_vertices, unsigned nb_faces, const unsigned int* labels,
    const float* vertices2, const unsigned int* indices2, unsigned nb_vertices2, unsigned nb_faces2, const unsigned int* labels2,
    const char* frame_json, char* movements_json, char** tooth_strs, int is_upper, int debug, int have_teeth, int have_movement)
{
    FakeGumGenConfig config;
    config.vertices = vertices;
    config.indices = indices;
    config.labels = labels;
    config.nb_vertices = nb_vertices;
    config.nb_faces = nb_faces;
    config.upper = true;
    config.mt = false;
    config.have_teeth = have_teeth == 1;
    config.debug = debug == 1;
    config.frame_json = frame_json;

    FakeGumGenConfig config2;
    config2.vertices = vertices2;
    config2.indices = indices2;
    config2.labels = labels2;
    config2.nb_vertices = nb_vertices2;
    config2.nb_faces = nb_faces2;
    config2.upper = false;
    config2.mt = false;
    config2.have_teeth = have_teeth == 1;
    config2.debug = debug == 1;
    config2.frame_json = frame_json;

    using hrc = std::chrono::high_resolution_clock;
    auto start = hrc::now();
    gumgen = std::make_unique<FakeGumGen>(config);
    gumgen2 = std::make_unique<FakeGumGen>(config2);
    std::cout << "Init: " << std::chrono::duration_cast<std::chrono::milliseconds>(hrc::now() - start).count() << "ms" << std::endl;
    return 0;
}

int EMSCRIPTEN_KEEPALIVE CreateStep( 
    float* transforms, 
    float** out_vertices, float** out_normals, int* out_vertex_sizes,
    float** out_vertices2, float** out_normals2, int* out_vertex_sizes2 )
{
    using hrc = std::chrono::high_resolution_clock;
    auto start = hrc::now();

    std::array<Eigen::Matrix4f, 49> mats;
    for(int i = 0; i < 49; i++)
        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                mats[i](k, j) = transforms[i * 16 + j * 4 + k];
    //std::fill(mats.begin(), mats.end(), Eigen::Matrix4f::Identity());

    Polyhedron mesh = gumgen->GenerateOneStep(mats);

    size_t nb_faces = mesh.size_of_facets();
    out_vertices[0] = new float[nb_faces * 3 * 3];
    out_normals[0] = new float[nb_faces * 3 * 3];
    out_vertex_sizes[0] = nb_faces * 3;

    int v_count = 0;
    for(auto hf = mesh.facets_begin(); hf != mesh.facets_end(); hf++)
    {
        Polyhedron::Traits::Kernel::Point_3 p[3] = {
                hf->halfedge()->vertex()->point(),
                hf->halfedge()->next()->vertex()->point(),
                hf->halfedge()->prev()->vertex()->point()
        };
        Polyhedron::Traits::Kernel::Vector_3 n[3] = {
                n[0] = hf->halfedge()->cgal_n,
                n[1] = hf->halfedge()->next()->cgal_n,
                n[2] = hf->halfedge()->prev()->cgal_n
        };

        for(int j = 0; j < 3; j++)
        {
            out_vertices[0][v_count * 3 + 0] = p[j].x();
            out_vertices[0][v_count * 3 + 1] = p[j].y();
            out_vertices[0][v_count * 3 + 2] = p[j].z();
            out_normals[0][v_count * 3 + 0] = n[j].x();
            out_normals[0][v_count * 3 + 1] = n[j].y();
            out_normals[0][v_count * 3 + 2] = n[j].z();
            v_count++;
        }
    }

    Polyhedron mesh2 = gumgen2->GenerateOneStep(mats);

    size_t nb_faces2 = mesh2.size_of_facets();
    out_vertices2[0] = new float[nb_faces2 * 3 * 3];
    out_normals2[0] = new float[nb_faces2 * 3 * 3];
    out_vertex_sizes2[0] = nb_faces2 * 3;

    int v_count2 = 0;
    for(auto hf = mesh2.facets_begin(); hf != mesh2.facets_end(); hf++)
    {
        Polyhedron::Traits::Kernel::Point_3 p[3] = {
                hf->halfedge()->vertex()->point(),
                hf->halfedge()->next()->vertex()->point(),
                hf->halfedge()->prev()->vertex()->point()
        };
        Polyhedron::Traits::Kernel::Vector_3 n[3] = {
                n[0] = hf->halfedge()->cgal_n,
                n[1] = hf->halfedge()->next()->cgal_n,
                n[2] = hf->halfedge()->prev()->cgal_n
        };

        for(int j = 0; j < 3; j++)
        {
            out_vertices2[0][v_count2 * 3 + 0] = p[j].x();
            out_vertices2[0][v_count2 * 3 + 1] = p[j].y();
            out_vertices2[0][v_count2 * 3 + 2] = p[j].z();
            out_normals2[0][v_count2 * 3 + 0] = n[j].x();
            out_normals2[0][v_count2 * 3 + 1] = n[j].y();
            out_normals2[0][v_count2 * 3 + 2] = n[j].z();
            v_count2++;
        }
    }
    std::cout << "Step : " << std::chrono::duration_cast<std::chrono::milliseconds>(hrc::now() - start).count() << "ms" << std::endl;
    return 0;
}


void Run( std::string uppermesh, std::string upperlabel, std::string lowermesh, std::string lowerlabel, std::string frame, 
std::string outuppermesh, std::string outlowermesh, std::string outupperweights, std::string outlowerweights)
{
    std::ifstream labels_ifs(upperlabel);
    std::string labels_str{ std::istreambuf_iterator<char>(labels_ifs), std::istreambuf_iterator<char>() };
    std::ifstream labels_ifs2(lowerlabel);
    std::string labels_str2{ std::istreambuf_iterator<char>(labels_ifs2), std::istreambuf_iterator<char>() };
    std::ifstream frame_ifs(frame);
    std::string frame_str{std::istreambuf_iterator<char>(frame_ifs), std::istreambuf_iterator<char>()};

    FakeGumGenConfig config;
    config.path_scanmesh = uppermesh;
    config.str_labels = labels_str;
    config.upper = true;
    config.mt = false;
    config.have_teeth = false;
    config.debug = false;
    config.frame_json = frame_str.c_str();

    FakeGumGenConfig config2;
    config2.path_scanmesh = lowermesh;
    config2.str_labels = labels_str2;
    config2.upper = false;
    config2.mt = false;
    config2.have_teeth = false;
    config2.debug = false;
    config2.frame_json = frame_str.c_str();

    gumgen = std::make_unique<FakeGumGen>(config);
    gumgen2 = std::make_unique<FakeGumGen>(config2);

    std::array<Eigen::Matrix4f, 49> mats;
    std::fill(mats.begin(), mats.end(), Eigen::Matrix4f::Identity());
    Polyhedron mesh = gumgen->GenerateOneStep(mats);
    Polyhedron mesh2 = gumgen2->GenerateOneStep(mats);
    
    std::ofstream ofs(outuppermesh, std::ios::binary);
    CGAL::IO::set_binary_mode(ofs);
    CGAL::IO::write_PLY(ofs, mesh);
    std::ofstream ofs2(outlowermesh, std::ios::binary);
    CGAL::IO::set_binary_mode(ofs2);
    CGAL::IO::write_PLY(ofs2, mesh2);
    gumgen->WriteWeights(outupperweights);
    gumgen2->WriteWeights(outlowerweights);
}

#ifndef __EMSCRIPTEN__
int main(int argc, char* argv[])
{

    // Run( "../../test/2/oral_scan_U.ply", "../../test/2/oral_scan_U.json", "../../test/2/oral_scan_L.ply", "../../test/2/oral_scan_L.json", "../../test/2/crown_frame.json",
    // "../../test/2/gumU.ply", "../../test/2/gumL.ply", "../../test/2/weightsU.json", "../../test/2/weightsL.json" );
    try
    {
        Run(std::string(argv[1]), std::string(argv[2]),std::string(argv[3]), std::string(argv[4]),
        std::string(argv[5]), std::string(argv[6]), std::string(argv[7]), std::string(argv[8]), std::string(argv[9]));
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }

    return 0;
}
#endif
}

#if 0
#include <pybind11/pybind11.h>
using py = pybind11;
PYBIND11_MODULE(gumgen, m)
{
    m.def("run", &Run, "run gumgen"),
    py::arg("uppermesh"), py::arg("upperlabel"),
    py::arg("lowermesh"), py::arg("lowerlabel"),
    py::arg("frame"),
    py::arg("outuppermesh"),
    py::arg("outlowermesh"),
    py::arg("outupperweights"),
    py::arg("outlowerweights")
}
#endif
