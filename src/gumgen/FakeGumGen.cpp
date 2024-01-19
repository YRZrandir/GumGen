#define NOMINMAX
#include "FakeGumGen.h"
#include <algorithm>
#include <chrono>
#include <thread>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/subdivision_method_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/io.h>
#include <Eigen/Geometry>
#include <Eigen/src/Geometry/Quaternion.h>
#include "../util/MathTypeConverter.h"
#include "../util/MeshFix.h"
#include "../eigen_bezier/eigen_bezier.hpp"
#include <igl/copyleft/cgal/mesh_boolean.h>
#ifndef __EMSCRIPTEN__
#include "tinycolormap.hpp"
#endif
class Timer
{
public:
    Timer() {}
    void start() { _start_time = std::chrono::high_resolution_clock::now(); }
    void stop() { _end_time = std::chrono::high_resolution_clock::now(); }
    int time() { return static_cast<int>(std::chrono::duration_cast<std::chrono::milliseconds>(_end_time - _start_time).count()); }

    std::chrono::time_point<std::chrono::high_resolution_clock> _start_time;
    std::chrono::time_point<std::chrono::high_resolution_clock> _end_time;
};

constexpr static bool EXPORT_RESAMPLE_GUM_LINES = false;
constexpr static bool EXPORT_SINGLE_MESH    = false;
constexpr static bool EXPORT_GUM_LINES      = false;
constexpr static bool EXPORT_BOTTOM_LINES   = false;
constexpr static bool EXPORT_SCANMESH       = false;
constexpr static bool EXPORT_CURVES         = false;

void write_gumline(const std::vector<Eigen::Vector3f>& polylines, std::string path, bool cycle = true)
{
#ifndef __EMSCRIPTEN__
    std::ofstream ofs1{std::string(path)};
    for(int gi = 0; gi < polylines.size(); gi++)
    {
        auto& p = polylines[gi];
        float c = (float)(gi) / polylines.size();
        auto color = tinycolormap::GetJetColor(c);
        
        ofs1 << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << ' ' << 0 << ' ' << 0 << ' ' << 0 << '\n';
    }
    int end = static_cast<int>(polylines.size());
    if(!cycle)
        end -= 1;
    for(int gi = 0; gi < end; gi++)
    {
        ofs1 << "l " << gi + 1 << ' ' << (gi + 2 > polylines.size() ? 1 : gi + 2) << '\n';
    }
    ofs1 << '\n';
    ofs1.close();
#endif
};


FakeGumGen::FakeGumGen( FakeGumGenConfig config )
    : _config( config )
{
    _scanmesh = std::make_unique<OralScanMesh>( config.upper );
#ifdef __EMSCRIPTEN__
    _scanmesh->LoadFromData(config.vertices, config.indices, config.nb_vertices, config.nb_faces);
    _scanmesh->ReadLabels(config.labels);
    _scanmesh->VertexLabelToFaceLabel();
#else
    if(!CGAL::IO::read_polygon_mesh(_config.path_scanmesh, *_scanmesh))
    {
        throw std::runtime_error("Cannot read input mesh!");
    }
    //_scanmesh->LoadFromString(_config.str_scanmesh);
    _scanmesh->ReadJsonLabels(_config.str_labels);
    _scanmesh->VertexLabelToFaceLabel();
#endif
    if(config.debug)
        std::cout << "Scanmesh: v=" << _scanmesh->size_of_vertices() << std::endl;
    OralScanMesh ori_scanmesh = *_scanmesh;

    _config.gum_width_list[11] = _config.gum_width_list[21] = _config.gum_width_list[31] = _config.gum_width_list[41] = 3.5f;
    _config.gum_width_list[12] = _config.gum_width_list[22] = _config.gum_width_list[32] = _config.gum_width_list[42] = 4.0f;
    _config.gum_width_list[13] = _config.gum_width_list[23] = _config.gum_width_list[33] = _config.gum_width_list[43] = 4.0f;
    _config.gum_width_list[14] = _config.gum_width_list[24] = _config.gum_width_list[34] = _config.gum_width_list[44] = 4.0f;
    _config.gum_width_list[15] = _config.gum_width_list[25] = _config.gum_width_list[35] = _config.gum_width_list[45] = 4.5f;
    _config.gum_width_list[16] = _config.gum_width_list[26] = _config.gum_width_list[36] = _config.gum_width_list[46] = 4.5f;
    _config.gum_width_list[17] = _config.gum_width_list[27] = _config.gum_width_list[37] = _config.gum_width_list[47] = 5.f;
    _config.gum_width_list[18] = _config.gum_width_list[28] = _config.gum_width_list[38] = _config.gum_width_list[48] = 5.f;

    _config.gum_depth_list[11] = _config.gum_depth_list[21] = _config.gum_depth_list[31] = _config.gum_depth_list[41] = 0.8f;
    _config.gum_depth_list[12] = _config.gum_depth_list[22] = _config.gum_depth_list[32] = _config.gum_depth_list[42] = 0.9f;
    _config.gum_depth_list[13] = _config.gum_depth_list[23] = _config.gum_depth_list[33] = _config.gum_depth_list[43] = 1.0f;
    _config.gum_depth_list[14] = _config.gum_depth_list[24] = _config.gum_depth_list[34] = _config.gum_depth_list[44] = 1.3f;
    _config.gum_depth_list[15] = _config.gum_depth_list[25] = _config.gum_depth_list[35] = _config.gum_depth_list[45] = 1.4f;
    _config.gum_depth_list[16] = _config.gum_depth_list[26] = _config.gum_depth_list[36] = _config.gum_depth_list[46] = 1.6f;
    _config.gum_depth_list[17] = _config.gum_depth_list[27] = _config.gum_depth_list[37] = _config.gum_depth_list[47] = 2.0f;
    _config.gum_depth_list[18] = _config.gum_depth_list[28] = _config.gum_depth_list[38] = _config.gum_depth_list[48] = 2.0f;

    Timer timer;
    timer.start();
    std::cout << "Preprocessing scan mesh..." << std::endl;
    _scanmesh->RemoveGumPart();
    //_scanmesh->ComputeOrientation();
    //_scanmesh->AlignMesh();
    if(EXPORT_SCANMESH)
    {
        CGAL::IO::write_PLY("scanmesh_toothpart" + std::to_string(_config.upper) + ".ply",  *_scanmesh);
    }
    auto aabb = _scanmesh->ComputeAABB();
    _config.gum_bottom = static_cast<float>(_config.upper ? (aabb.zmax() + 10.f) : (aabb.zmin() - 10.f));
    timer.stop();
    std::cout << "finish. "<< timer.time() << std::endl;

    if(_config.debug)
    {
        std::cout << "Create Teeth Meshes...";
        timer.start();
    }
    _teeth_meshes = std::make_unique<TeethMeshes>(_config.upper);
    _teeth_meshes->InitFromOralScanMesh(*_scanmesh);
    //_teeth_meshes->AlignMeshes(_scanmesh->_obb_up, _scanmesh->_obb_center);
    _teeth_meshes->LoadTeethOritFromJson(_config.frame_json);
    if(_config.debug)
    {
        timer.stop();
        std::cout << "Finish..." << timer.time() << std::endl;
    }

    if(_config.debug)
    {
        std::cout << "Query Edges...";
        timer.start();
    }
    std::array<std::vector<hHalfedge>, 49> gumline_edges = GetGumlineEdges();
    if(_config.debug)
        std::cout << "Finish. " << timer.time() << std::endl;
    
    if(_config.debug)
    {
        std::cout << "Sorting Edges...";
        timer.start();
    }
    for (int i = 10; i < gumline_edges.size(); ++i)
    {
        if(!gumline_edges[i].empty())
        {
            std::cout << "Label " << i << ": " << gumline_edges[i].size() << ";";
            SortGumlineEdges(gumline_edges[i]);
        }
    }
    if(_config.debug)
    {
        timer.stop();
        std::cout << "Finish." << timer.time() << std::endl;
    }

    if(_config.debug)
    {
        timer.start();
        std::cout << "Resample and Smooth Edges...";
    }
    
    AABBTree aabb_tree(CGAL::faces(ori_scanmesh).first, CGAL::faces(ori_scanmesh).second, ori_scanmesh);
    if(aabb_tree.empty())
    {
        std::cout << "error" << ori_scanmesh.size_of_vertices() << " " << ori_scanmesh.size_of_facets() << std::endl;
    }
    for (int i = 11; i < 49; i++)
    {
        if (gumline_edges[i].empty())
            continue;
        for (auto hh : gumline_edges[i])
            _gumlines[i].push_back( ToEigen( hh->vertex()->point() ).cast<float>() );
        if(EXPORT_GUM_LINES)
            write_gumline(_gumlines[i], std::string("gumline_old") + std::to_string(i) + std::string(".obj"));
        _gumlines[i] = ResampleGumLine(_gumlines[i]);
        //SmoothAndProj(_gumlines[i], 10, aabb_tree);
        auto center = std::accumulate(_gumlines[i].begin(), _gumlines[i].end(), Eigen::Vector3f(0.f, 0.f, 0.f)) / _gumlines[i].size();
        //_gumlines[i] = ResampleGumLine( _gumlines[i], center, _teeth_meshes->_teeth_pca_dir[i] );
    }
    if(_config.debug)
    {
        timer.stop();
        std::cout << "Finish." << timer.time() << std::endl;
    }

    std::array<std::optional<Eigen::Vector3f>, 49> centers = ToothCenterPoints( _gumlines );

    //FillEmptyTeeth(_gumlines, centers);
    if(EXPORT_RESAMPLE_GUM_LINES)
    {
        for(int i = 11; i < 49; i++)
        {
            if(!_gumlines[i].empty())
            {
                write_gumline(_gumlines[i], std::string("gumline") + std::to_string(i) + std::string(".obj"));
            }
        }
    }
    if(_config.upper)
    {
        for (int i = 18; i >= 11; i--)
            if (centers[i].has_value())
                _center_pairs.emplace_back( centers[i].value(), i  );
        for (int i = 21; i <= 27; i++)
            if (centers[i].has_value())
                _center_pairs.emplace_back( centers[i].value(), i  );

    }
    else
    {
        for (int i = 37; i >= 31; i--)
            if (centers[i].has_value())
                _center_pairs.emplace_back( centers[i].value(), i  );
        for (int i = 41; i <= 47; i++)
            if (centers[i].has_value())
                _center_pairs.emplace_back( centers[i].value(), i  );
    }

    if(_config.debug)
    {
        std::cout << "Creating gum for " << _center_pairs.size() << " teeth" << std::endl;
    }

}

void
FakeGumGen::AlignGumMeshBackward(Polyhedron& mesh )
{
    Eigen::Quaternion<double> q = Eigen::Quaternion<double>::FromTwoVectors( _scanmesh->_obb_up, Eigen::Vector3d( 0, 0, _config.upper ? -1 : 1 ) );
    Eigen::Matrix3d rot = q.matrix();
    rot.transposeInPlace();

    for (auto hv = mesh.vertices_begin(); hv != mesh.vertices_end(); hv++)
    {
        auto p = ToEigen( hv->point() );
        p -= _scanmesh->_obb_center;
        p = rot * p;
        p += _scanmesh->_obb_center;
        hv->point() = CGAL::Point_3<CGAL::Epick>( p.x(), p.y(), p.z() );
    }

    for(auto hh = mesh.halfedges_begin(); hh != mesh.halfedges_end(); hh++)
    {
        auto n = ToEigen( hh->cgal_n );
        n = rot * n;
        hh->cgal_n = ToCGAL<double, Polyhedron::Traits::Kernel>(n);
    }
}

Polyhedron
FakeGumGen::GenerateOneStep( const std::array<Eigen::Matrix4f, 49>& transforms )
{
    std::array<std::vector<Eigen::Vector3f>, 49> cur_gumline = _gumlines;
    std::vector<std::pair<Eigen::Vector3f, int>> cur_center_pairs = _center_pairs;
    std::array<Eigen::Vector3f, 49> cur_pca_dirs = _teeth_meshes->_teeth_pca_dir;
    std::array<Eigen::Vector3f, 49> cur_pca_centroids = _teeth_meshes->_teeth_pca_centroid;
    std::vector<std::unique_ptr<Polyhedron>> meshes(cur_center_pairs.size());

    //Transform meshes according to movement path info.
    for (int l = 0; l < 49; l++)
    {
        auto& M = transforms[l];
        for (auto& p : cur_gumline[l])
        {
            Eigen::Vector4f pw(p.x(), p.y(), p.z(), 1.0);
            pw = M * pw;
            p = Eigen::Vector3f(pw.x(), pw.y(), pw.z());
        }


        Eigen::Vector4f wc = M * Eigen::Vector4f(cur_pca_centroids[l].x(), cur_pca_centroids[l].y(), cur_pca_centroids[l].z(), 1.0);
        cur_pca_centroids[l] = Eigen::Vector3f(wc.x(), wc.y(), wc.z());
        Eigen::Vector4f wd = M * Eigen::Vector4f(cur_pca_dirs[l].x(), cur_pca_dirs[l].y(), cur_pca_dirs[l].z(), 0.0);
        cur_pca_dirs[l] = Eigen::Vector3f(wd.x(), wd.y(), wd.z());
    }

    for (auto& pair : cur_center_pairs)
    {
        Eigen::Vector4f wct = transforms[pair.second] * Eigen::Vector4f(pair.first.x(), pair.first.y(), pair.first.z(), 1.0);
        pair.first = Eigen::Vector3f(wct.x(), wct.y(), wct.z());
    }

    float new_bottom_pos =  _config.gum_bottom;

    using Vertex_point_map = boost::property_map<Polyhedron, CGAL::vertex_point_t>::type;
    using TreeTraitsBase = CGAL::Search_traits_3<Polyhedron::Traits::Kernel>;
    using TreeTraits = CGAL::Search_traits_adapter<boost::graph_traits<Polyhedron>::vertex_descriptor, Vertex_point_map, TreeTraitsBase>;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree KdTree;

    std::vector<Eigen::MatrixXd> igl_Vs(cur_center_pairs.size());
    std::vector<Eigen::MatrixXi> igl_Fs(cur_center_pairs.size());
    std::vector<std::pair<std::unique_ptr<Polyhedron>, int>> seperate_gums(cur_center_pairs.size());
    std::vector<std::pair<std::unique_ptr<KdTree>, Vertex_point_map>> seperate_gum_acc(cur_center_pairs.size());
#pragma omp parallel for
    for (int i = 0; i < cur_center_pairs.size(); i++)
    {
        int label = cur_center_pairs[i].second;
        Eigen::Vector3f dir;
        std::optional<Eigen::Vector3f> to_prev_center;
        std::optional<Eigen::Vector3f> to_next_center;
        if (i == 0)
        {
            dir = (cur_center_pairs[1].first - cur_center_pairs[0].first).normalized();
            to_next_center = cur_center_pairs[1].first - cur_center_pairs[0].first;
        }
        else if (i == cur_center_pairs.size() - 1)
        {
            dir = (cur_center_pairs[i].first - cur_center_pairs[i - 1].first).normalized();
            to_prev_center = cur_center_pairs[i - 1].first - cur_center_pairs[i].first;
        }
        else
        {
            dir = (cur_center_pairs[i].first - cur_center_pairs[i - 1].first).normalized();
            dir += (cur_center_pairs[i + 1].first - cur_center_pairs[i].first).normalized();
            dir /= 2.f;
            to_next_center = cur_center_pairs[i + 1].first - cur_center_pairs[i].first;
            to_prev_center = cur_center_pairs[i - 1].first - cur_center_pairs[i].first;
        }

        auto gum = CreateGumForOneTooth( cur_gumline[label], cur_center_pairs[i].first, dir, label,
         cur_pca_centroids[label], cur_pca_dirs[label], new_bottom_pos, to_prev_center, to_next_center );
        if(EXPORT_SINGLE_MESH)
        {
            gum->WriteOBJ(std::string("gum") + std::to_string(label) + std::string(".obj"));
        }
        if(_config.debug)
        {
            std::cout << "Validate Gum...";
            std::cout << "is_closed()=";
            std::cout << gum->is_closed();

            std::cout << ". is_valid()=";
            std::cout << gum->is_valid();

            std::cout << ". degenerate faces=";
            std::vector<hFacet> degenerate_faces;
            CGAL::Polygon_mesh_processing::degenerate_faces( *gum, std::back_inserter(degenerate_faces));
            std::cout << degenerate_faces.size() << std::endl;

            std::cout << "non-manifold vertices=";
            std::vector<hHalfedge> non_manifold_vertices;
            CGAL::Polygon_mesh_processing::non_manifold_vertices( *gum, std::back_inserter(non_manifold_vertices));
            std::cout << non_manifold_vertices.size() << std::endl;
            printf("Single Gum %d, v=%zd\n", label, gum->size_of_vertices());
        }

        auto [V, F] = gum->WriteToEigen();
        igl_Vs[i] = std::move(V);
        igl_Fs[i] = std::move(F);
        Vertex_point_map vpmap = CGAL::get(CGAL::vertex_point, *gum);
        seperate_gum_acc[i] = std::make_pair(std::make_unique<KdTree>(CGAL::vertices(*gum).begin(), CGAL::vertices(*gum).end(), KdTree::Splitter(), TreeTraits(vpmap)), vpmap);
        seperate_gums[i] = std::move(std::make_pair(std::move(gum), label));
    }

    Polyhedron final_gum;
    Eigen::MatrixXd final_V(0, 3);
    Eigen::MatrixXi final_F(0, 3);
    for(int i = 0; i < igl_Vs.size(); i++)
    {
        if(i == 0)
        {
            final_V = std::move(igl_Vs[i]);
            final_F = std::move(igl_Fs[i]);
        }
        else
        {
            std::cout << "try union..." << std::endl;
            Eigen::MatrixXd tempV;
            Eigen::MatrixXi tempF;
            igl::copyleft::cgal::mesh_boolean(igl_Vs[i], igl_Fs[i], final_V, final_F, igl::MESH_BOOLEAN_TYPE_UNION, tempV, tempF);
            final_V = std::move(tempV);
            final_F = std::move(tempF);

            // fix non-manifold
            std::vector<Polyhedron::Point_3> final_points;
            std::vector<TTriangle<Polyhedron::Vertex::size_type>> final_faces;
            for(int r = 0; r < final_V.rows(); r++)
                final_points.emplace_back(final_V(r, 0), final_V(r, 1), final_V(r, 2));
            for(int r = 0; r < final_F.rows(); r++)
                final_faces.emplace_back(final_F(r, 0), final_F(r, 1), final_F(r, 2));
            Polyhedron fixed_mesh;
            FixMesh(final_points, final_faces, fixed_mesh, true, 1000, false, false, 0, 0, false, 10);
            auto pair = fixed_mesh.WriteToEigen();
            final_V = std::move(pair.first);
            final_F = std::move(pair.second);
        }
    }
    final_gum.LoadFromEigen(final_V, final_F);

    std::vector<hHalfedge> border_edges = final_gum.FindHoles();

    for (auto hh : border_edges)
    {
        final_gum.CloseHole(hh, true, false);
    }
    if(_config.debug)
    {
        std::cout << "Complete. " << std::endl;
    }
    if(_config.debug)
    {
        std::cout << "Validate final mesh...";
        std::cout << "is_closed()=";
        std::cout << final_gum.is_closed();

        std::cout << ". is_valid()=";
        std::cout << final_gum.is_valid();

        std::cout << ". degenerate faces=";
        std::vector<hFacet> degenerate_faces;
        CGAL::Polygon_mesh_processing::degenerate_faces( final_gum, std::back_inserter(degenerate_faces));
        std::cout << degenerate_faces.size() << std::endl;

        std::cout << "non-manifold vertices=";
        std::vector<hHalfedge> non_manifold_vertices;
        CGAL::Polygon_mesh_processing::non_manifold_vertices( final_gum, std::back_inserter(non_manifold_vertices));
        std::cout << non_manifold_vertices.size() << std::endl;
    }

    // Mark labels
    for(auto& hv : CGAL::vertices(final_gum))
    {
        double min_dist = std::numeric_limits<double>::max();
        for(auto& [gum, vpmap] : seperate_gum_acc)
        {
            Neighbor_search::Distance dist(vpmap);
            Neighbor_search search(*gum, hv->point(), 1, 0.0, true, dist);
            if(search.begin()->second < min_dist)
            {
                min_dist = search.begin()->second;
                hv->label = search.begin()->first->label;
            }
        }
    }
    final_gum.WriteOBJ("labeled_gum" + std::to_string(_config.upper) + ".obj");

    std::vector<hFacet> smooth_faces;
    for(auto hf = final_gum.facets_begin(); hf != final_gum.facets_end(); hf++)
    {
        if(hf->halfedge()->vertex()->label != 1 && hf->halfedge()->prev()->vertex()->label != 1 && hf->halfedge()->next()->vertex()->label != 1)
        {
            smooth_faces.push_back(hf);
        }
    }
    CGAL::Polygon_mesh_processing::angle_and_area_smoothing( smooth_faces, final_gum, CGAL::Polygon_mesh_processing::parameters::number_of_iterations( 1 ).use_safety_constraints( false ) );
    //CGAL::Subdivision_method_3::Loop_subdivision( *_final_gums[step], CGAL::parameters::number_of_iterations( 1 ) );
    //CGAL::Polygon_mesh_processing::tangential_relaxation(CGAL::vertices(*_final_gums[step]), *_final_gums[step], CGAL::Polygon_mesh_processing::parameters::number_of_iterations(5));
    LaplacianSmooth( final_gum, 5, true );

    ClipFinalGum( final_gum, new_bottom_pos );
    FixMesh(final_gum, true, 1000, false, true, 0, 0.f, false, 10);
    if(_config.debug)
    {
        std::cout << "Validate final mesh...";
        std::cout << "is_closed()=";
        std::cout << final_gum.is_closed();

        std::cout << ". is_valid()=";
        std::cout << final_gum.is_valid();

        std::cout << ". degenerate faces=";
        std::vector<hFacet> degenerate_faces;
        CGAL::Polygon_mesh_processing::degenerate_faces( final_gum, std::back_inserter(degenerate_faces));
        std::cout << degenerate_faces.size() << std::endl;

        std::cout << "non-manifold vertices=";
        std::vector<hHalfedge> non_manifold_vertices;
        CGAL::Polygon_mesh_processing::non_manifold_vertices( final_gum, std::back_inserter(non_manifold_vertices));
        std::cout << non_manifold_vertices.size() << std::endl;
    }

    _weights = ComputeWeights( final_gum );

    //AlignGumMeshBackward( final_gum );
    final_gum.UpdateNormal();
    return final_gum;
}

void
FakeGumGen::ClipFinalGum( Polyhedron& mesh, float gum_bottom)
{
    using Kernel = Polyhedron::Traits::Kernel;
    Kernel::Plane_3 clip_plane{
        Kernel::Point_3{0.0, 0.0, gum_bottom + (_config.upper ? -3.0 : 3.0)},
        Kernel::Vector_3{0.0, 0.0, (_config.upper ? 1 : -1)}
    };

    CGAL::Polygon_mesh_processing::clip(mesh, clip_plane, CGAL::Polygon_mesh_processing::parameters::allow_self_intersections(true));
    auto border_edges = mesh.FindHoles();
    for(auto hh : border_edges)
    {
        auto [patch_faces, patch_vertices] = mesh.CloseHole(hh, false, false);
//        for(auto hf : patch_faces)
//        {
//            auto p0 = ToEigen(hf->halfedge()->vertex()->point());
//            auto p1 = ToEigen(hf->halfedge()->next()->vertex()->point());
//            auto p2 = ToEigen(hf->halfedge()->prev()->vertex()->point());
//            Eigen::Vector3d n = ((p1 - p0).cross(p2 - p0));
//            n = Eigen::Vector3d(0.0, 0.0, _config.upper ? 1.0 : -1.0);
//            hf->halfedge()->cgal_n = ToCGAL<double, Polyhedron::Traits::Kernel>(n);
//            hf->halfedge()->next()->cgal_n = ToCGAL<double, Polyhedron::Traits::Kernel>(n);
//            hf->halfedge()->prev()->cgal_n = ToCGAL<double, Polyhedron::Traits::Kernel>(n);
//        }
    }
}

std::array<std::vector<hHalfedge>, 49>
FakeGumGen::GetGumlineEdges()
{
    std::array<std::vector<hHalfedge>, 49> gumline_edges;
    for (auto hh = _scanmesh->halfedges_begin(); hh != _scanmesh->halfedges_end(); hh++)
    {
        if (hh->opposite() == nullptr)
            continue;
        if (hh->is_border() && !hh->opposite()->is_border() && hh->opposite()->facet()->label != 0)
        {
            gumline_edges[hh->opposite()->facet()->label ].push_back( hh );
        }
    }

    for(int i = 11; i < 49; i++)
        if(gumline_edges[i].size() != 0)
            std::cout << "GumLine of label " << i << " = " << gumline_edges[i].size() << std::endl;

    return gumline_edges;
}

void
FakeGumGen::SortGumlineEdges(std::vector<hHalfedge>& tooth_edges )
{
    std::vector<std::list<hHalfedge>> sort_edges( 1 );
    sort_edges[0].push_back( tooth_edges.front() );
    std::unordered_set<hHalfedge> all_edges( std::next(tooth_edges.begin()), tooth_edges.end() );

    while (!all_edges.empty())
    {
        hHalfedge last_edge = sort_edges.back().back();
        hHalfedge front_edge = sort_edges.back().front();
        bool found = false;
        for (auto it = all_edges.begin(); it != all_edges.end(); it++)
        {
            if (last_edge->vertex() == (*it)->prev()->vertex())
            {
                sort_edges.back().push_back( *it );
                all_edges.erase( it );
                found = true;
                break;
            }
            else if (front_edge->prev()->vertex() == (*it)->vertex())
            {
                sort_edges.back().push_front( *it );
                all_edges.erase( it );
                found = true;
                break;
            }
        }
        if (!found)
        {
            sort_edges.push_back( std::list<hHalfedge>() );
            sort_edges.back().push_back( *all_edges.begin() );
            all_edges.erase( all_edges.begin() );
        }
    }
    if (sort_edges.size() > 2)
    {
        sort_edges.erase( std::remove_if( sort_edges.begin(), sort_edges.end(), []( std::list<hHalfedge>& edges ) { return edges.size() <= 3; } ), sort_edges.end() );
    }
    if (sort_edges.size() == 1)
    {
        if(_config.debug)
            std::cout << "1 segment." << std::endl;
        tooth_edges = std::vector<hHalfedge>( sort_edges[0].begin(), sort_edges[0].end() );
    }
    else if (sort_edges.size() == 2)
    {
        if(_config.debug)
            std::cout << "2 segments." << std::endl;
        tooth_edges.clear();
        for (const auto& li : sort_edges)
            tooth_edges.insert( tooth_edges.end(), li.begin(), li.end() );
    }
    else
    {
        if(_config.debug)
            std::cout << sort_edges.size() << "segments." << std::endl;
        Polyhedron::Traits::Kernel::Vector_3 center(0.f, 0.f, 0.f);
        double w_sum = 0.f;
        for(hHalfedge hh : tooth_edges)
        {
            double w = std::sqrt((hh->vertex()->point() - hh->prev()->vertex()->point()).squared_length()) * 0.5;
            w_sum += w;
            center += (hh->vertex()->point() - CGAL::ORIGIN) * w;
            center += (hh->prev()->vertex()->point() - CGAL::ORIGIN) * w;
        }
        center /= w_sum * 2.f;
        Eigen::Vector3d ptcenter = ToEigen(center);

        tooth_edges.clear();
        tooth_edges.insert( tooth_edges.end(), sort_edges[0].begin(), sort_edges[0].end() );
        sort_edges.erase( sort_edges.begin() );

        while (!sort_edges.empty())
        {
            Eigen::Vector3d p_last = ToEigen(tooth_edges.back()->vertex()->point());
            Eigen::Vector3d d_last = (p_last - ptcenter).normalized();
            Eigen::Vector3d vec_last = p_last - ToEigen(tooth_edges.back()->prev()->vertex()->point());

            auto min_it = std::max_element( sort_edges.begin(), sort_edges.end(),
                [&]( std::list<hHalfedge>& lh, std::list<hHalfedge>& rh )
                {
                    Eigen::Vector3d p1 = ToEigen(lh.front()->vertex()->point());
                    Eigen::Vector3d p2 = ToEigen(rh.front()->vertex()->point());
                    
                    Eigen::Vector3d vec1 = (ToEigen(lh.front()->vertex()->point()) - p_last);
                    Eigen::Vector3d vec2 = (ToEigen(rh.front()->vertex()->point()) - p_last);

                    double dot1 = vec1.normalized().dot(vec_last);
                    double dot2 = vec2.normalized().dot(vec_last);

                    return dot1 < dot2;
                } );
            tooth_edges.insert( tooth_edges.end(), min_it->begin(), min_it->end() );
            sort_edges.erase( min_it );
        }
    }
}

std::array<std::optional<Eigen::Vector3f>, 49>
FakeGumGen::ToothCenterPoints(const std::array<std::vector<Eigen::Vector3f>, 49>& points )
{
    std::array<std::optional<Eigen::Vector3f>, 49> centers;
    for (int i = 11; i < 49; i++)
    {
        if (!points[i].empty())
            centers[i] = std::accumulate( points[i].begin(), points[i].end(), Eigen::Vector3f( 0.f, 0.f, 0.f ) ) / static_cast<float>(points[i].size());
    }
    return centers;
}

void
FakeGumGen::FillEmptyTeeth(std::array<std::vector<Eigen::Vector3f>, 49>& gumlines, std::array<std::optional<Eigen::Vector3f>, 49>& centers)
{
    const std::vector<int> indices = _config.upper ?
     std::vector<int>{ 17, 16, 15, 14, 13, 12, 11, 21, 22, 23, 24, 25, 26, 27 } : 
     std::vector<int>{ 37, 36, 35, 34, 33, 32, 31, 41, 42, 43, 44, 45, 46, 47 };

    int start_idx = 0;
    while(!centers[indices[start_idx]].has_value())
        start_idx++;
    for(int idx = start_idx; idx < indices.size(); idx++)
    {
        int label = indices[idx];
        if(!centers[label].has_value())
        {
            
            int next_valid_idx = idx + 1;
            while(next_valid_idx < indices.size() && !centers[indices[next_valid_idx]].has_value())
                next_valid_idx++;
            if(next_valid_idx >= indices.size())
                break;

            auto prev_valid_pos = centers[indices[idx - 1]].value();
            auto next_valid_pos = centers[indices[next_valid_idx]].value();
            int diff = next_valid_idx - idx;
            Eigen::Vector3f pos = 1.f / (float)(diff + 1) * ( next_valid_pos - prev_valid_pos) + prev_valid_pos;
            centers[label] = pos;
            Eigen::Vector3f up = (_teeth_meshes->_teeth_pca_dir[indices[idx - 1]] + _teeth_meshes->_teeth_pca_dir[indices[next_valid_idx]]) * 0.5f;
            _teeth_meshes->_teeth_pca_dir[indices[idx]] = up;
            _teeth_meshes->_teeth_updir[indices[idx]] = up;
            _teeth_meshes->_teeth_pca_centroid[indices[idx]] = (_teeth_meshes->_teeth_pca_centroid[indices[idx - 1]] + _teeth_meshes->_teeth_pca_centroid[indices[next_valid_idx]]) * 0.5f;

            for(int j = 0; j < _config.nb_sample; j++)
            {
                float angle = (float)(_config.upper ? j : (_config.nb_sample - j)) / _config.nb_sample * 2.f * 3.1415f;
                Eigen::Vector3f u = up.cross(Eigen::Vector3f(1, 0, 0)).normalized();
                Eigen::Vector3f v = up.cross(u).normalized();
                gumlines[label].push_back(centers[label].value() + _config.gum_width_list[label] * 0.5f * (std::cos(angle) * v + std::sin(angle) * u));
            }

            std::cout << "Add Tooth Between " << indices[idx - 1] << " and " << indices[next_valid_idx] << std::endl;
        }
    }
}

std::vector<Eigen::Vector3f>
FakeGumGen::ResampleGumLine(const std::vector<Eigen::Vector3f>& points )
{
    float gumline_len = 0.f;
    for (int i = 0; i < points.size(); i++)
    {
        gumline_len += (points[i] - points[(i + 1) % points.size()]).norm();
    }
    float step_len = gumline_len / _config.nb_sample;

    std::vector<Eigen::Vector3f> gumline_resample;
    gumline_resample.push_back( points[0] );
    int iedge = 0;
    float cur_len = 0.f;
    for (int i = 0; i < _config.nb_sample; i++)
    {
        //Find next sample
        float sampleline_len = gumline_resample.size() * step_len;
        float nextsample_len = sampleline_len + step_len;
        while (nextsample_len > cur_len)
        {
            iedge++;
            if (iedge >= points.size())
                break;

            cur_len += (points[iedge] - points[(iedge + 1) % points.size()]).norm();
        }
        if (iedge >= points.size())
            break;
        float len_diff = cur_len - nextsample_len;
        Eigen::Vector3f p0 = points[iedge];
        Eigen::Vector3f p1 = points[(iedge + 1) % points.size()];
        Eigen::Vector3f p_sample = p1 + (p0 - p1).normalized() * len_diff;
        gumline_resample.push_back( p_sample );
    }

    return gumline_resample;
}

std::vector<Eigen::Vector3f>
FakeGumGen::ResampleGumLine( const std::vector<Eigen::Vector3f>& edges, Eigen::Vector3f center, Eigen::Vector3f updir )
{
    Eigen::Vector3f _u = updir.cross(Eigen::Vector3f(1, 0, 0)).normalized();
    Eigen::Vector3f _v = updir.cross(_u).normalized();
    KernelEpick::Vector_3 c(center.x(), center.y(), center.z());
    KernelEpick::Vector_3 up(updir.x(), updir.y(), updir.z());

    std::vector<Eigen::Vector3f> results;
    for(int i = 0; i < _config.nb_sample; i++)
    {
        float ang = (1.0f - (float)(i) / _config.nb_sample) * 2 * 3.14159f;
        float sin = std::sin(ang);
        float cos = std::cos(ang);
        
        Eigen::Vector3f plane_d = _u * cos + _v * sin;
        Eigen::Vector3f plane_p = center + plane_d * 1000.0;
        KernelEpick::Triangle_3 tri(CGAL::ORIGIN + ToCGAL<float, KernelEpick>(plane_p) , CGAL::ORIGIN + c - up * 1000.0, CGAL::ORIGIN + c + up * 1000.0);
        std::vector<Eigen::Vector3f> intersect_points;
        
        for(int j = 0; j < edges.size(); j++)
        {
            int idx0 = j;
            int idx1 = (j + 1) % edges.size();
            auto p0 = edges[idx0];
            auto p1 = edges[idx1];
            
            auto ret = CGAL::intersection(KernelEpick::Segment_3(CGAL::ORIGIN + ToCGAL<float, KernelEpick>(p0), CGAL::ORIGIN + ToCGAL<float, KernelEpick>(p1)), tri);
            if(ret.has_value() && ret->which() == 0)
            {
                auto p = boost::get<KernelEpick::Point_3>(*ret);
                Eigen::Vector3f pn = (ToEigen(p).cast<float>() - center).dot(updir) * updir;
                Eigen::Vector3f pt = (ToEigen(p).cast<float>() - center) - pn;
                intersect_points.push_back(ToEigen(p).cast<float>());
            }
        }
        if(!intersect_points.empty())
        {
            Eigen::Vector3f pos = *std::max_element(intersect_points.begin(), intersect_points.end(), [&center, &updir](Eigen::Vector3f& left, Eigen::Vector3f& right) {
                float dist0 = (left - center).dot(updir);
                float dist1 = (right - center).dot(updir);
                return dist0 < dist1;
            });
            results.push_back(pos);
        }
    }

    return results;
}

void
FakeGumGen::SmoothGumLine(std::vector<Eigen::Vector3f>& points, int nb_iteration )
{
    for (int ite = 0; ite < nb_iteration; ite++)
    {
        std::vector<Eigen::Vector3f> result( points.size() );
        for (int i = 0; i < points.size(); i++)
        {
            size_t prev = (i == 0) ? points.size() - 1 : i - 1;
            size_t next = (i == points.size() - 1) ? 0 : i + 1;
            result[i] = (points[prev] + points[next]) * 0.5f;
        }
        points = std::move( result );
    }
}

void
FakeGumGen::SmoothAndProj( std::vector<Eigen::Vector3f>& points, int nb_iterate, const FakeGumGen::AABBTree& aabb_tree )
{
    for(int i = 0; i < nb_iterate; i++)
    {
        std::vector<Eigen::Vector3f> result( points.size() );
        for(int j = 0; j < points.size(); j++)
        {
            size_t prev = (j == 0) ? points.size() - 1 : j - 1;
            size_t next = (j == points.size() - 1) ? 0 : j + 1;
            result[j] = (points[prev] + points[next]) * 0.5f;
            result[j] = ToEigen(aabb_tree.closest_point(CGAL::ORIGIN + ToCGAL<float, KernelEpick>(result[j]))).cast<float>();
        }
        points = std::move(result);
    }
}

std::unique_ptr<Polyhedron>
FakeGumGen::CreateGumForOneTooth(const std::vector<Eigen::Vector3f>& points, Eigen::Vector3f center,
 Eigen::Vector3f dir, int label, Eigen::Vector3f pca_centroid, Eigen::Vector3f pca_dir,
  float gum_bottom, std::optional<Eigen::Vector3f> to_prev_center, std::optional<Eigen::Vector3f> to_next_center )
{
    if(_config.debug)
    {
        std::cout << "Try generate gum for tooth: " << label << std::endl;
    }
    std::unique_ptr<Polyhedron> gum = std::make_unique<Polyhedron>( false );

    //dir.z() = 0.f;
    dir.normalize();

    Eigen::Vector3f center_bottom = center + pca_dir * (_config.upper ? (center.z() - gum_bottom) : (gum_bottom - center.z()));
    Eigen::Vector3f side_dir = (pca_dir.cross( dir )).normalized();

    /*
    7  0  1
     \ | /
    6- c -2
     / | \
    5  4  3
    */

    std::array<Eigen::Vector2d, 16> outline_points2d = {
        Eigen::Vector2d{0.0, 1.0},
        Eigen::Vector2d{0.3, 1.0},
        Eigen::Vector2d{1.0, 1.0},
        Eigen::Vector2d{1.0, 0.3},
        Eigen::Vector2d{1.0, 0.0},
        Eigen::Vector2d{1.0, -0.3},
        Eigen::Vector2d{1.0, -1.0},
        Eigen::Vector2d{0.3, -1.0},
        Eigen::Vector2d{0.0, -1.0},
        Eigen::Vector2d{-0.3, -1.0},
        Eigen::Vector2d{-1.0, -1.0},
        Eigen::Vector2d{-1.0, -0.3},
        Eigen::Vector2d{-1.0, 0.0},
        Eigen::Vector2d{-1.0, 0.3},
        Eigen::Vector2d{-1.0, 1.0},
        Eigen::Vector2d{-0.3, 1.0},
    };
    for(auto& p2d : outline_points2d)
    {
        if(p2d.y() > 0)
            p2d.y() *= _config.gum_width_list[label] * 1.5f;
        else
            p2d.y() *= _config.gum_width_list[label] * 0.8f;
        //p2d.x() *= _config.gum_width_list[label] * 0.7f;
        if(p2d.x() > 0)
        {
            if(to_next_center.has_value())
                p2d.x() *= to_next_center.value().norm() * 2.0;
        }
        else
        {
            if(to_prev_center.has_value())
                p2d.x() *= to_prev_center.value().norm() * 2.0;
        }
    }

    std::vector<std::shared_ptr<Bezier::Curve>> outline_curve_list;
    for (int i = 1; i < outline_points2d.size(); i += 2)
    {
        Eigen::MatrixX2d points_2d( 3, 2 );
        points_2d.row( 0 ) = outline_points2d[i];
        points_2d.row( 1 ) = outline_points2d[(i + 1) % outline_points2d.size()];
        points_2d.row( 2 ) = outline_points2d[(i + 2) % outline_points2d.size()];
        outline_curve_list.push_back( std::make_shared<Bezier::Curve>( points_2d ) );
    }
    Bezier::PolyCurve outline_curve( outline_curve_list );

    std::vector<Eigen::Vector3f> outline_resample( points.size() );
    Bezier::PointVector polyline = outline_curve.polyline( 1.0001, 0.1f );
    CGAL::Simple_cartesian<float>::Plane_3 bottom_proj_plane(CGAL::ORIGIN + ToCGAL(center_bottom), ToCGAL(pca_dir));
    Eigen::Vector3f u = dir.normalized();
    Eigen::Vector3f v = side_dir.normalized();

    for (size_t i = 0; i < points.size(); i++)
    {
        Eigen::Vector3f proj_diff = ToEigen(bottom_proj_plane.projection(CGAL::ORIGIN + ToCGAL(points[i]))) - center_bottom;
        Eigen::Vector2f gumline2d( proj_diff.dot(u), proj_diff.dot(v) );
        Eigen::Vector2f center2d( 0.f, 0.f );
        for (int j = 0; j < polyline.size(); j++)
        {
            Eigen::Vector2f seg0 = polyline[j].cast<float>();
            Eigen::Vector2f seg1 = polyline[(j + 1) % polyline.size()].cast<float>();

            Eigen::Vector2f o1 = center2d;
            Eigen::Vector2f d1 = (gumline2d - center2d).normalized();
            Eigen::Vector2f o2 = seg0;
            Eigen::Vector2f d2 = (seg1 - seg0).normalized();

            Eigen::Matrix2f m;
            m( 0, 0 ) = d1.x();
            m( 0, 1 ) = -d2.x();
            m( 1, 0 ) = d1.y();
            m( 1, 1 ) = -d2.y();
            if (std::abs( m.determinant() ) < 1e-5f)
                continue;
            Eigen::Vector2f y = o2 - o1;
            Eigen::Vector2f x = m.inverse() * y;
            float t = x[0];
            float s = x[1];

            if (t > 0 && s > 0 && s < (seg1 - seg0).norm())
            {
                Eigen::Vector2f p2d = o2 + s * d2;
                outline_resample[i] = p2d.x() * u + p2d.y() * v + center_bottom;
                break;
            }
        }
       // outline_resample[i] = polyline[min_j].x() * u + polyline[min_j].y() * v + center_bottom;
    }
    if(EXPORT_BOTTOM_LINES)
        write_gumline(outline_resample, std::string("bottomline") + std::to_string(label) + std::string(".obj"));

    //Gen side face
    Eigen::Vector3f center_ori = center;
    center -= static_cast<float>(_config.gum_depth_list[label]) * 3 * pca_dir;
    std::vector<std::vector<Eigen::Vector3f>> side_points( points.size() );
    std::vector<int> side_points_labels( points.size() );
    for (int i = 0; i < points.size(); i++)
    {
        Eigen::Vector3f gumline_point = points[i];
        Eigen::Vector3f center_to_gumline = gumline_point - center;
        Eigen::Vector3f center_to_gumline_v = pca_dir.dot(center_to_gumline) * pca_dir;
        Eigen::Vector3f center_to_gumline_h = center_to_gumline - center_to_gumline_v;

        //gumline_point.z() += (_config.upper ? -1.f : 1.f) * 0.5f;
        //float thickness = 1.0f * std::abs(center_to_gumline_h.dot(dir));
        float thickness = 1.0f;
        // if(to_next_center.has_value())
        // {
        //     if(center_to_gumline_h.dot(to_next_center.value()) > 0)
        //         thickness += std::pow(center_to_gumline_h.normalized().dot(to_next_center.value().normalized()), 2.0) * to_next_center.value().norm();
        // }
        // else
        // {
        //     thickness += std::pow(center_to_gumline_h.normalized().dot(to_prev_center.value().normalized()), 2.0) * to_prev_center.value().norm();
        // }
        // if(to_prev_center.has_value())
        // {
        //     if(center_to_gumline_h.dot(to_prev_center.value()) > 0)
        //         thickness += std::pow(center_to_gumline_h.normalized().dot(to_prev_center.value().normalized()), 2.0) * to_prev_center.value().norm();
        // }
        // else
        // {
        //     thickness += std::pow(center_to_gumline_h.normalized().dot(to_next_center.value().normalized()), 2.0) * to_next_center.value().norm();
        // }

        Eigen::Vector3f outline_point = outline_resample[i];
        Eigen::Vector3f ctrl_cg = gumline_point - center_to_gumline_h * 0.5f;
        Eigen::Vector3f ctrl_cg1 = center + center_to_gumline_h * 0.5f;
        Eigen::Vector3f extrude_dir = center_to_gumline_h.normalized();
        //extrude_dir = ((gumline_point - center).cross(_teeth_updir[label])).cross(_teeth_updir[label]).normalized();
        float extrude_w = (5.f + _config.gum_width_list[label] * 0.5f) * 0.1f * _config.extrude_width;
        float h = gumline_point.z() - outline_point.z();
        Eigen::Vector3f ctrl_go0 = outline_point + (center - center_bottom) + center_to_gumline_h.normalized() * thickness;
        Eigen::Vector3f ctrl_go1 = outline_point + center_to_gumline_h.normalized() * thickness;
        //ctrl_go0 += extrude_dir * extrude_w * 1.0f;
        //ctrl_go1 += extrude_dir * extrude_w * 1.0f;
        
        eigen_bezier::Curve<float> curve1({center, ctrl_cg1, ctrl_cg, gumline_point});
        eigen_bezier::Curve<float> curve2({gumline_point, ctrl_go0, ctrl_go1, outline_point});
        for(int j = 0; j < 10; j++)
        {
            float t = (float)j / 10.f;
            side_points[i].push_back(curve1.Eval(t));
        }
        for(int j = 0; j < 10; j++)
        {
            float t = (float)j / 10.f;
            // if(j < 5)
            //     t = j / 20.f;
            // else
            //     t = 0.25f + (j - 5.f) / (20.f / 3);
            side_points[i].push_back(curve2.Eval(t));
        }
        side_points_labels[i] = label;

        if(EXPORT_CURVES)
        {
            write_gumline(side_points[i], "side" + std::to_string(label) + ".obj", false);
        }
    }

    if(_config.debug)
    {
        printf("Try Deletgate...\n");
    }
    FakeGumBuilder2<Polyhedron::HalfedgeDS> builder( side_points );
    gum->delegate( builder );
    if(_config.debug)
    {
        printf("Cleaning...\n");
    }
    std::vector<std::pair<hFacet, hFacet>> intersected_tris;
    CGAL::Polygon_mesh_processing::self_intersections( CGAL::faces( *gum ), *gum, std::back_inserter( intersected_tris ) );
    std::unordered_set<hFacet> erase_facets_set;
    for (const auto& pair : intersected_tris)
    {
        erase_facets_set.insert( pair.first );
        erase_facets_set.insert( pair.second );
    }
    for (auto hf : erase_facets_set)
    {
        gum->erase_facet( hf->halfedge() );
    }
    erase_facets_set.clear();
    gum->keep_largest_connected_components( 1 );

    // std::vector<hHalfedge> non_manifold_vertices;
    // CGAL::Polygon_mesh_processing::non_manifold_vertices( *gum, std::back_inserter(non_manifold_vertices));
    // for(auto hh : non_manifold_vertices)
    // {
    //     for(auto hf : CGAL::faces_around_target(hh, *gum))
    //     {
    //         erase_facets_set.insert(hf);
    //     }
    // }
    // for (auto hf : erase_facets_set)
    // {
    //     gum->erase_facet( hf->halfedge() );
    // }
    
    FixMeshWithLabel(*gum, false, 1000, false, true, 0, 0.f, false, 10);

    if(_config.debug)
    {
        printf("close holes...");
    }
    std::vector<hHalfedge> border_edges;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles( *gum, std::back_inserter( border_edges ) );
    for (auto h : border_edges)
    {
        std::vector<hFacet> patch_facets;
        std::vector<hVertex> patch_vertices;
        CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(*gum, h, std::back_inserter(patch_facets), std::back_inserter(patch_vertices));
    }
    //LaplacianSmooth(*gum, 5, true);
    //gum->WriteOFF("./test" + std::to_string(label) + ".off");

    FixMeshWithLabel(*gum, true, 1000, false, false, 0, 0.f, false, 10);
    if(_config.debug)
    {
        std::cout << "Success" << std::endl;
    }
    
    return gum;
}

void
FakeGumGen::LaplacianSmooth(Polyhedron& mesh, int iterate_num, bool skip_gumline )
{
    for (int i = 0; i < iterate_num; i++)
    {
        for (auto hv = mesh.vertices_begin(); hv != mesh.vertices_end(); hv++)
        {
            if (hv->degree() == 0)
                continue;
            if (skip_gumline && hv->label == 1)
                continue;

            float x = 0.f;
            float y = 0.f;
            float z = 0.f;

            float w_sum = 0.f;
            for(auto hv2 : CGAL::vertices_around_target(hv, mesh))
            {
                float w = 1.f;
                x += static_cast<float>(hv2->point().x()) * w;
                y += static_cast<float>(hv2->point().y()) * w;
                z += static_cast<float>(hv2->point().z()) * w;
                w_sum += w;
            }

            x /= w_sum;
            y /= w_sum;
            z /= w_sum;

            if (w_sum < 1e-7f)
            {
                continue;
            }

            hv->point() = { x, y, z };
        }
    }
}

std::tuple<Eigen::Vector3f, Eigen::Vector3f, Eigen::Vector3f>
FakeGumGen::ComputeUpdirAndCenterPointsFromOBB(const std::array<Eigen::Vector3f, 8>& obb )
{
    Eigen::Vector3f candidate0 = obb[1] - obb[0];
    Eigen::Vector3f candidate1 = obb[3] - obb[0];
    Eigen::Vector3f candidate2 = obb[5] - obb[0];
    Eigen::Vector3f test{0, 0, (_config.upper ? -1.f : 1.f)};
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

std::vector<std::array<std::pair<int, float>, 3>>
FakeGumGen::ComputeWeights( const Polyhedron& gum_mesh )
{
    std::vector<std::array<std::pair<int, float>, 3>> weight_table;
    for(auto hv = gum_mesh.vertices_begin(); hv != gum_mesh.vertices_end(); hv++)
    {
        Eigen::Vector3f p = ToEigen(hv->point()).cast<float>();
        
        std::vector<std::pair<int, float>> weights;
        for(int t = 11; t < 49; t++)
        {
            if(_gumlines[t].empty())
                continue;
            auto tcenter = _teeth_meshes->_teeth_pca_centroid[t];
            auto tdir = _teeth_meshes->_teeth_pca_dir[t];
            auto& borders = _gumlines[t];

            //min dist to borders
            double min_dist_to_borders = std::numeric_limits<float>::max();
            for(int i = 0; i < borders.size(); i++)
            {
                int i0 = i;
                int i1 = (i + 1) % borders.size();
                auto e0 = _gumlines[t][i0];
                auto e1 = _gumlines[t][i1];
                double dist_to_segment = CGAL::squared_distance(
                    KernelEpick::Segment_3{ CGAL::ORIGIN + ToCGAL<float, KernelEpick>(e0), CGAL::ORIGIN + ToCGAL<float, KernelEpick>(e1) },
                    CGAL::ORIGIN + ToCGAL<float, KernelEpick>(p));
                if(dist_to_segment < min_dist_to_borders)
                    min_dist_to_borders;
            }

            double dist_to_frame = CGAL::squared_distance(
                KernelEpick::Line_3{ CGAL::ORIGIN + ToCGAL<float, KernelEpick>(tcenter), ToCGAL<float, KernelEpick>(tdir) },
                CGAL::ORIGIN + ToCGAL<float, KernelEpick>(p)
            );

            double final_dist = std::min(min_dist_to_borders, dist_to_frame);
            weights.emplace_back(t, static_cast<float>(final_dist));
        }
        
        std::sort(weights.begin(), weights.end(), []( auto& lh, auto& rh ) { return lh.second < rh.second; });
        std::array<std::pair<int, float>, 3> weight_table_this;
        weight_table_this.fill({0, 0});
        if(weights[0].second < 1e-8f)
        {
            weight_table_this[0] = {weights[0].first, 1.0};
        }
        else
        {
            float w_sum = 0.f;
            for(int i = 0; i < 3 && i < weights.size(); i++)
            {
                weight_table_this[i] = {weights[i].first, 1.f / weights[i].second};
                w_sum += weight_table_this[i].second;
            }
            for(auto& [label, w] : weight_table_this)
                if(label != 0)
                    w /= w_sum;
        }
        weight_table.push_back(weight_table_this);
    }
    std::cout << "Weights: " << weight_table.size() << std::endl;
    return weight_table;
}

void FakeGumGen::WriteWeights( std::string path )
{
    using namespace nlohmann;
    json js;
    js["weights"] = _weights;
    std::ofstream ofs( path );
    ofs << js;
    //ofs.write(reinterpret_cast<const char*>(_weights.data()), sizeof(_weights[0]) * _weights.size());
    //ofs.close();
}