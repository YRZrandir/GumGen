//
// Created by yrz on 3/19/23.
//
#define NOMINMAX
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <unordered_set>
#include <thread>
#include "gumgen_support_classes.h"
#include "nlohmann/json.hpp"
#include "../util/MathTypeConverter.h"
#include "../util/MeshFix.h"

OralScanMesh::OralScanMesh(bool upper)
    :Polyhedron(false), _upper(upper)
{

}

void OralScanMesh::ReadJsonLabels( const std::string& str )
{
    using namespace nlohmann;
    std::stringstream label_ifs( str );
    json data = json::parse( label_ifs );
    if (data.find( "labels" ) == data.end())
    {
        std::cout << "Invalid Json" << std::endl;
        return;
    }
    std::vector<int> labels = data["labels"].get<std::vector<int>>();
    int count = 0;
    for (auto hv = vertices_begin(); hv != vertices_end(); hv++)
    {
        hv->label = labels[count++];
        if (hv->label == 100 || hv->label % 10 == 8)
            hv->label = 0;
    }
}

void OralScanMesh::ReadLabels( const unsigned int* labels)
{
    int count = 0;
    for(auto hv = vertices_begin(); hv != vertices_end(); hv++)
    {
        hv->label = labels[count++];
        if(hv->label == 100 || hv->label % 10 == 8)
            hv->label = 0;
    }

    std::unordered_map<unsigned int, int> counts;
    for(auto hv = vertices_begin(); hv != vertices_end(); hv++)
    {
        counts[hv->label]++;
    }

    for(auto& [label, num] : counts)
    {
        std::cout << "Label " << label << " = " << num << std::endl;
    }
}

void OralScanMesh::VertexLabelToFaceLabel()
{
    for(auto hf = facets_begin(); hf != facets_end(); hf++)
    {
        int l0 = hf->halfedge()->vertex()->label;
        int l1 = hf->halfedge()->next()->vertex()->label;
        int l2 = hf->halfedge()->prev()->vertex()->label;

        if(l0 == l1 && l0 == l2)
            hf->label = l0;
        else
        {
            if(l0 >= l1 && l0 >= l2)
                hf->label = l0;
            else if(l1 >= l2 && l1 >= l0)
                hf->label = l1;
            else
                hf->label = l2;
        }
    }
}

void OralScanMesh::ComputeOrientation()
{
    auto obb_points = ComputeOBB();

    std::array<Eigen::Vector3d, 8> obb_points2;
    for (int i = 0; i < 8; i++)
    {
        obb_points2[i] = ToEigen( obb_points[i] );
        _obb_center += obb_points2[i];
    }
    _obb_center /= 8.0;

    Eigen::Vector3d candidate0 = ((obb_points2[1] - obb_points2[0]).normalized().cross( obb_points2[2] - obb_points2[0] ).normalized());
    Eigen::Vector3d candidate1 = ((obb_points2[4] - obb_points2[0]).normalized().cross( obb_points2[5] - obb_points2[0] ).normalized());
    Eigen::Vector3d candidate2 = ((obb_points2[1] - obb_points2[0]).normalized().cross( obb_points2[6] - obb_points2[0] ).normalized());
    Eigen::Vector3d test{0, 0, (_upper ? -1.0 : 1.0)};
    double dot0 = candidate0.dot(test);
    double dot1 = candidate1.dot(test);
    double dot2 = candidate2.dot(test);
    Eigen::Vector3d up{};
    if(std::abs(dot0) > std::abs(dot1) && std::abs(dot0) > std::abs(dot2))
    {
        if(dot0 < 0)
            up = -candidate0;
        else
            up = candidate0;
    }
    else if(std::abs(dot1) > std::abs(dot2))
    {
        if(dot1 < 0)
            up =  -candidate1;
        else
            up = candidate1;
    }
    else
    {
        if(dot2 < 0)
            up =  -candidate2;
        else
            up =  candidate2;
    }
    _obb_up = up;

    auto aabb = ComputeAABB();
    _obb_center = Eigen::Vector3d(aabb.xmax() + aabb.xmin(), aabb.ymax() + aabb.ymin(), aabb.zmax() + aabb.zmin());
    _obb_center *= 0.5;

    std::cout << "ScanMesh: " << _obb_up.x() << ',' << _obb_up.y() << ',' << _obb_up.z() << "; " << _obb_center.x() << ',' << _obb_center.y() << ',' << _obb_center.z() << std::endl;
}

void OralScanMesh::AlignMesh()
{
    Eigen::Quaternion<double> q = Eigen::Quaternion<double>::FromTwoVectors( _obb_up, Eigen::Vector3d( 0, 0, _upper ? -1 : 1 ) );
    Eigen::Matrix3d rot = q.matrix();
    for (auto hv = vertices_begin(); hv != vertices_end(); hv++)
    {
        auto p = ToEigen( hv->point() );
        p -= _obb_center;
        p = rot * p;
        p += _obb_center;
        hv->point() = CGAL::Point_3<CGAL::Epick>( p.x(), p.y(), p.z() );
    }
}

void OralScanMesh::AlignMeshInverse()
{
    Eigen::Quaternion<double> q = Eigen::Quaternion<double>::FromTwoVectors( _obb_up, Eigen::Vector3d( 0, 0, _upper ? -1 : 1 ) );
    Eigen::Matrix3d rot = q.matrix();
    rot.transposeInPlace();
    for (auto hv = vertices_begin(); hv != vertices_end(); hv++)
    {
        auto p = ToEigen( hv->point() );
        p -= _obb_center;
        p = rot * p;
        p += _obb_center;
        hv->point() = CGAL::Point_3<CGAL::Epick>( p.x(), p.y(), p.z() );
    }
}

void OralScanMesh::RemoveGumPart()
{
    std::vector<hHalfedge> facets_to_remove;
    for (auto hf = facets_begin(); hf != facets_end(); hf++)
    {
        if(hf->label == 0)
        {
            facets_to_remove.push_back( hf->halfedge() );
        }
    }
    for (auto hh : facets_to_remove)
    {
        erase_facet( hh );
    }

    //CGAL::Polygon_mesh_processing::keep_largest_connected_components( *this, 1 );

    auto hole_edges = FindHoles();
    std::vector<std::vector<hHalfedge>> holes;
    for(hHalfedge hh : hole_edges)
        holes.push_back(HoleEdges(hh));

    std::sort( holes.begin(), holes.end(), []( auto& left, auto& right ) { return left.size() < right.size(); } );
    holes.pop_back();
    for (auto& hole : holes)
    {
        if(hole.size() <= 100)
            CloseHole(hole[0], false);
    }
}

TeethMeshes::TeethMeshes(bool upper)
    :_upper(upper)
{
}

void TeethMeshes::InitFromOralScanMesh(OralScanMesh& oral_scan_mesh)
{
    //Partition the oral scan mesh
    CGAL::Face_filtered_graph<OralScanMesh> filtered_mesh(oral_scan_mesh);
    std::array<std::vector<hFacet>, 49> faces_group_by_labels;
    for(auto hf = oral_scan_mesh.facets_begin(); hf != oral_scan_mesh.facets_end(); hf++)
    {
        int l0 = hf->halfedge()->vertex()->label;
        int l1 = hf->halfedge()->next()->vertex()->label;
        int l2 = hf->halfedge()->prev()->vertex()->label;
        int facelabel = l0;
        if(l0 == l1)
            facelabel = l0;
        else if(l1 == l2 && l1 != l0)
            facelabel = l1;
        else if(l1 != l2 && l1 != l0)
            facelabel = 0;
        if( facelabel != 0 )
        {
            faces_group_by_labels[hf->halfedge()->vertex()->label].push_back(hFacet(hf));
        }
    }

    for(int i = 11; i < 49; i++)
    {
        std::vector<hFacet>& selected_faces = faces_group_by_labels[i];
        if(!selected_faces.empty())
        {
            _teeth[i] = std::make_unique<Polyhedron>(false);
            filtered_mesh.set_selected_faces( selected_faces );
            CGAL::copy_face_graph(filtered_mesh, *_teeth[i]);
            FixMesh(*_teeth[i], true, 1000, false, true, 0, 0.f, false, 10);
            CGAL::Polygon_mesh_processing::keep_largest_connected_components( *_teeth[i], 1 );
        }
    }
    // for(int i = 11; i < 49; i++)
    // {
    //     if(_teeth[i] == nullptr)
    //         continue;

    //     auto edges = _teeth[i]->FindHoles();
    //     for(auto hh : edges)
    //     {
    //         auto[new_faces, new_vertices] = _teeth[i]->CloseHole(hh, true, false);
    //         std::unordered_set<hFacet> smooth_faces{new_faces.begin(), new_faces.end()};
    //         std::unordered_set<hFacet> additional_faces;
    //         for(hFacet hf : smooth_faces)
    //         {
    //             for(hFacet f : CGAL::faces_around_face(hf->halfedge(), *_teeth[i]))
    //             {
    //                 if(f != nullptr && smooth_faces.count(f) == 0)
    //                     additional_faces.insert(f);
    //             }
    //         }
    //         for(hFacet hf : additional_faces)
    //         {
    //             for(hFacet f : CGAL::faces_around_face(hf->halfedge(), *_teeth[i]))
    //             {
    //                 if(f != nullptr && smooth_faces.count(f) == 0)
    //                     additional_faces.insert(f);
    //             }
    //         }
    //         for(hFacet hf : additional_faces)
    //         {
    //             for(hFacet f : CGAL::faces_around_face(hf->halfedge(), *_teeth[i]))
    //             {
    //                 if(f != nullptr && smooth_faces.count(f) == 0)
    //                     additional_faces.insert(f);
    //             }
    //         }
    //         smooth_faces.insert(additional_faces.begin(), additional_faces.end());
    //         std::vector<hVertex> smooth_vertices;
    //         for(hFacet hf : smooth_faces)
    //         {
    //             if(hf == nullptr)
    //                 continue;
    //             smooth_vertices.push_back(hf->halfedge()->vertex());
    //             smooth_vertices.push_back(hf->halfedge()->next()->vertex());
    //             smooth_vertices.push_back(hf->halfedge()->prev()->vertex());
    //         }
    //         _teeth[i]->LaplacianSmooth(smooth_vertices, 5, 0);
    //         //CGAL::Polygon_mesh_processing::angle_and_area_smoothing(smooth_faces, *_teeth[i], CGAL::Polygon_mesh_processing::parameters::number_of_iterations(1) );
    //     }

        //_teeth[i]->WriteOFF("./tooth" + std::to_string(i) + ".off");
    //}
}

void TeethMeshes::InitFromFile(const std::array<std::string, 49>& paths)
{
    for(int i = 10; i < 49; i++)
    {
        if(i % 10 >= 1 && i % 10 <= 7 && !paths[i].empty())
        {
            _teeth[i] = std::make_unique<Polyhedron>(paths[i], false);
        }
    }
}

void TeethMeshes::InitFromString( const std::array<std::string, 49>& strings)
{
    for(int i = 10; i < 49; i++)
    {
        if(i % 10 >= 1 && i % 10 <= 7 && !strings[i].empty())
        {
            _teeth[i] = std::make_unique<Polyhedron>(strings[i]);
        }
    }
}

void TeethMeshes::AlignMeshes( Eigen::Vector3d obb_up, Eigen::Vector3d obb_center)
{
    _obb_up = obb_up;
    _obb_center = obb_center;

    Eigen::Quaternion<float> q = Eigen::Quaternion<float>::FromTwoVectors( _obb_up.cast<float>(), Eigen::Vector3f( 0.f, 0.f, _upper ? -1.f : 1.f ) );
    Eigen::Matrix3f rot = q.matrix();

    for(auto& d : _teeth_updir)
        d = rot * d;
    for(auto& d : _teeth_pca_dir)
        d = rot * d;
    for(auto& p : _teeth_pca_centroid)
        p = (rot * (p - _obb_center.cast<float>())) + _obb_center.cast<float>();

//#pragma omp parallel for default(none) shared(rot)
    for(auto& tooth : _teeth)
    {
        if(tooth != nullptr)
        {
            for (auto hv = tooth->vertices_begin(); hv != tooth->vertices_end(); hv++)
            {
                Eigen::Vector3f p = ToEigen( hv->point() ).cast<float>();
                p -= _obb_center.cast<float>();
                p = rot * p;
                p += _obb_center.cast<float>();
                hv->point() = CGAL::Point_3<CGAL::Epick>( p.x(), p.y(), p.z() );
            }
        }
    }

}

void TeethMeshes::AlignMeshesInverse()
{
    Eigen::Quaternion<double> q = Eigen::Quaternion<double>::FromTwoVectors( _obb_up, Eigen::Vector3d( 0, 0, _upper ? -1 : 1 ));
    Eigen::Matrix3d rot = q.matrix();
    rot.transposeInPlace();
    for(auto& tooth : _teeth)
    {
        if(tooth != nullptr)
        {
            for (auto hv = tooth->vertices_begin(); hv != tooth->vertices_end(); hv++)
            {
                auto p = ToEigen( hv->point() );
                p -= _obb_center;
                p = rot * p;
                p += _obb_center;
                hv->point() = CGAL::Point_3<CGAL::Epick>( p.x(), p.y(), p.z() );
            }
        }
    }
}

void TeethMeshes::ComputeTeethOBB()
{
    // std::vector<std::thread> threads;
    // for(size_t i = 0; i < _teeth.size(); i++)
    // {
    //     if(_teeth[i] != nullptr)
    //     {
    //         threads.emplace_back([&, this, i]() {
    //             auto obb_points_0 = _teeth[i]->ComputeOBB();
    //             std::array<Eigen::Vector3f, 8> obb_points;
    //             for(size_t j = 0; j < 8; j++)
    //                 obb_points[j] = ToEigen(obb_points_0[j]).cast<float>();
    //             _teeth_obb[i] = obb_points;
    //             auto [updir, top_center, bottom_center] = ComputeUpdirAndCenterPointsFromOBB(obb_points);
    //             _teeth_updir[i] = updir;
    //             _teeth_bottom_centers[i] = bottom_center;
    //         });
    //     }
    // }
    // for(auto& t : threads)
    //     t.join();
    for(size_t i = 0; i < _teeth.size(); i++)
    {
        if(_teeth[i] == nullptr)
            continue;

        auto obb_points_0 = _teeth[i]->ComputeOBB();
        std::array<Eigen::Vector3f, 8> obb_points;
        for(size_t j = 0; j < 8; j++)
            obb_points[j] = ToEigen(obb_points_0[j]).cast<float>();
        _teeth_obb[i] = obb_points;
        auto [updir, top_center, bottom_center] = ComputeUpdirAndCenterPointsFromOBB(obb_points);
        _teeth_updir[i] = updir;
        _teeth_bottom_centers[i] = bottom_center;
    }
}

void TeethMeshes::LoadTeethOritFromJson( const char* str )
{
    nlohmann::json data = nlohmann::json::parse(str);
    
    for(int i = 11; i < 49; i++)
    {
        if(i % 10 >= 8)
            continue;
        if(data.find(std::to_string(i)) != data.end())
        {
            std::vector<std::vector<float>> frame_raw = data[std::to_string(i)].get<std::vector<std::vector<float>>>();
            Eigen::Vector3f centroid{ frame_raw[3][0], frame_raw[3][1], frame_raw[3][2] };
            Eigen::Vector3f up{frame_raw[2][0], frame_raw[2][1], frame_raw[2][2]};
            _teeth_updir[i] = up;
            _teeth_pca_dir[i] = up;
            _teeth_pca_centroid[i] = centroid;
        }
    }
}

void TeethMeshes::ComputeTeethPCA()
{
    for(int i = 0; i < _teeth.size(); i++)
    {
        if(_teeth[i] == nullptr)
            continue;

        auto [line, centroid] = _teeth[i]->ComputePCA();
        _teeth_pca_dir[i] = ToEigen(line.direction()).cast<float>().normalized();
        _teeth_pca_centroid[i] = ToEigen(centroid).cast<float>();

        Eigen::Vector3f test{0.f, 0.f, (_upper ? -1.f : 1.f)};
        if(_teeth_pca_dir[i].dot(test) < 0)
            _teeth_pca_dir[i] = -_teeth_pca_dir[i];
    }
}

std::tuple<Eigen::Vector3f, Eigen::Vector3f, Eigen::Vector3f> TeethMeshes::ComputeUpdirAndCenterPointsFromOBB( const std::array<Eigen::Vector3f, 8>& obb)
{
    Eigen::Vector3f candidate0 = obb[1] - obb[0];
    Eigen::Vector3f candidate1 = obb[3] - obb[0];
    Eigen::Vector3f candidate2 = obb[5] - obb[0];
    Eigen::Vector3f test{0, 0, (_upper ? -1.f : 1.f)};
    float dot0 = std::abs(candidate0.dot(test));
    float dot1 = std::abs(candidate1.dot(test));
    float dot2 = std::abs(candidate2.dot(test));
    Eigen::Vector3f bottom1;
    Eigen::Vector3f bottom2;
    Eigen::Vector3f updir;
    Eigen::Vector3f bottom_center;
    if(dot0 > dot1 && dot0 > dot2)
    {
        bottom1 = (obb[0] + obb[4]) * 0.5f;
        bottom2 = (obb[1] + obb[7]) * 0.5f;
        updir = candidate0.normalized();
    }
    else if(dot1 > dot2)
    {
        bottom1 = (obb[0] + obb[6]) * 0.5f;
        bottom2 = (obb[3] + obb[7]) * 0.5f;
        updir = candidate1.normalized();
    }
    else
    {
        bottom1 = (obb[0] + obb[2]) * 0.5f;
        bottom2 = (obb[5] + obb[7]) * 0.5f;
        updir = candidate2.normalized();
    }
    if((bottom1 - bottom2).dot(test) > 0)
    {
        return std::make_tuple(updir, bottom1, bottom2);
    }
    else
    {
        return std::make_tuple(updir, bottom2, bottom1);
    }
}