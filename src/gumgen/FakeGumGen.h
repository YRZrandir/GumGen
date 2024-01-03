//
// Created by yrz on 12/21/22.
//

#ifndef ELASTICITY_FAKEGUMGEN_H
#define ELASTICITY_FAKEGUMGEN_H
#include "Polyhedron.h"
#include <fstream>
#include <nlohmann/json.hpp>
#include "../Bezier/bezier.h"
#include "../Bezier/polycurve.h"
#include "gumgen_support_classes.h"

struct FakeGumGenConfig
{
    std::string path_scanmesh;
    std::string path_labels;
    std::string str_labels;

    const float* vertices;
    const unsigned* indices;
    const unsigned* labels;
    unsigned nb_vertices;
    unsigned nb_faces;
    const char* frame_json;
    
    std::array<std::string, 49> str_teeth;
    bool debug = true;
    bool upper = false;
    bool mt = false;
    bool have_teeth = false;
    float gum_bottom = -45.f;
    std::array<float, 49> gum_width_list;
    std::array<float, 49> gum_depth_list;
    int nb_sample = 35;
    float extrude_width = 1.25f;
};

template <typename HDS>
class FakeGumBuilder2 : public CGAL::Modifier_base<HDS>
{
public:
    FakeGumBuilder2( const std::vector<std::vector<Eigen::Vector3f>>& points ) : _points( points ) {};
    virtual void operator()( HDS& hds ) override;

private:
    std::vector<std::vector<Eigen::Vector3f>> _points;
};

template <typename HDS>
void FakeGumBuilder2<HDS>::operator()( HDS& hds )
{
    CGAL::Polyhedron_incremental_builder_3<HDS> builder( hds, false );
    builder.begin_surface( _points.size() * _points[0].size(), _points.size() * (_points[0].size() - 1) * 2 - _points.size() );

    std::vector<std::vector<size_t>> indices;
    size_t count = 0;
    for (int i = 0; i < _points.size(); i++)
    {
        indices.push_back( std::vector<size_t>() );
        const auto& points = _points[i];
        for (int j = 0; j < points.size(); j++)
        {
            hVertex vh = builder.add_vertex( Point3( points[j].x(), points[j].y(), points[j].z() ) );
            if (j == points.size() / 2 || j == points.size() - 1)
                vh->label = 1;
            if (j == points.size() - 1)
                vh->label = 2;
            indices[i].push_back( count );
            count++;
        }
    }

    for (int i = 0; i < _points.size(); i++)
    {
        indices[i][0] = indices[0][0];
    }

    for (int i = 0; i < indices.size(); i++)
    {
        const std::vector<size_t>& points0 = indices[i];
        const std::vector<size_t>& points1 = indices[(i + 1) % _points.size()];
        for (int j = 0; j < points0.size() - 1; j++)
        {
            size_t left_up = points0[j];
            size_t left_bottom = points0[j + 1];
            size_t right_up = points1[j];
            size_t right_bottom = points1[j + 1];
            builder.begin_facet();
            builder.add_vertex_to_facet( left_up );
            builder.add_vertex_to_facet( right_bottom );
            builder.add_vertex_to_facet( left_bottom );
            builder.end_facet();
            if (j == 0)
                continue;
            builder.begin_facet();
            builder.add_vertex_to_facet( left_up );
            builder.add_vertex_to_facet( right_up );
            builder.add_vertex_to_facet( right_bottom );
            builder.end_facet();
        }
    }
    builder.remove_unconnected_vertices();
    builder.end_surface();
}


class FakeGumGen
{
public:
    explicit FakeGumGen(FakeGumGenConfig config );

    void AlignGumMeshBackward( Polyhedron& mesh );
    void AlignTeethMesh();
    std::tuple<Eigen::Vector3f, Eigen::Vector3f, Eigen::Vector3f> ComputeUpdirAndCenterPointsFromOBB( const std::array<Eigen::Vector3f, 8>& obb );
    Polyhedron GenerateOneStep( const std::array<Eigen::Matrix4f, 49>& transforms );
    void ClipFinalGum( Polyhedron& mesh, float gum_bottom );

    std::array<std::vector<hHalfedge>, 49> GetGumlineEdges();
    void SortGumlineEdges( std::vector<hHalfedge>& tooth_edges );
    std::vector<Eigen::Vector3f> ResampleGumLine( const std::vector<Eigen::Vector3f>& edges );
    std::vector<Eigen::Vector3f> ResampleGumLine( const std::vector<Eigen::Vector3f>& edges, Eigen::Vector3f center, Eigen::Vector3f updir );
    void SmoothGumLine( std::vector<Eigen::Vector3f>& points, int nb_iteration );
    std::array<std::optional<Eigen::Vector3f>, 49> ToothCenterPoints( const std::array<std::vector<Eigen::Vector3f>, 49>& edges );
    void FillEmptyTeeth( std::array<std::vector<Eigen::Vector3f>, 49>& gumlines, std::array<std::optional<Eigen::Vector3f>, 49>& centers);

    std::unique_ptr<Polyhedron> CreateGumForOneTooth(
        const std::vector<Eigen::Vector3f>& points, Eigen::Vector3f center, Eigen::Vector3f dir, int label,
         Eigen::Vector3f pca_centroid, Eigen::Vector3f pca_dir, float gum_bottom );

    void LaplacianSmooth( Polyhedron& mesh, int iterate_num, bool skip_gumline );
    std::vector<std::array<std::pair<int, float>, 3>> ComputeWeights( const Polyhedron& gum_mesh );
    void WriteWeights( std::string path );
    //std::vector<std::array<std::pair<int, float>, 3>> ReadWeights( void* weights );
    std::unique_ptr<OralScanMesh> _scanmesh;
    std::unique_ptr<TeethMeshes> _teeth_meshes;
    std::array<std::vector<Eigen::Vector3f>, 49> _gumlines;
    std::vector<std::pair<Eigen::Vector3f, int>> _center_pairs;

    std::array<std::vector<std::unique_ptr<Polyhedron>>, 22> _meshes;
    std::array<std::unique_ptr<Polyhedron>, 22> _final_gums;
    std::array<std::array<Eigen::Matrix4f, 49>, 22> _movements;
    std::array<Eigen::Vector3f, 49> _centroids;
    std::vector<std::array<std::pair<int, float>, 3>> _weights;

private:
    FakeGumGenConfig _config{};
};

#endif