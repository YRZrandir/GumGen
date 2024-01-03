//
// Created by yrz on 7/7/22.
//

#include <cmath>
#include <sstream>
#include "Polyhedron.h"
#include "../util/MathTypeConverter.h"
#include <nlohmann/json.hpp>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/linear_least_squares_fitting_3.h>

Polyhedron::Polyhedron( const std::string& path, bool render )
    :_render( render )
{

    
//    {
//        const static auto pt2vec3 = []( const tinyobj::real_t* p ) { return glm::vec3( p[0], p[1], p[2] );  };
//        const static auto pt2vec2 = []( const tinyobj::real_t* p ) { return glm::vec2( p[0], p[1] ); };
//        auto& materials = reader.GetMaterials();
//        MaterialInfos main_info;
//        if (!materials.empty())
//        {
//            const auto& main_mat = materials[0];
//            main_info.alpha = main_mat.transmittance[0];
//            main_info.diffuseColor = pt2vec3( main_mat.diffuse );
//            main_info.specularColor = pt2vec3( main_mat.specular );
//            main_info.diffuseTexPath = main_mat.diffuse_texname;
//            main_info.specularTexPath = main_mat.specular_texname;
//            main_info.normalTexPath = main_mat.bump_texname;
//            main_info.roughnessTexPath = main_mat.roughness_texname;
//            main_info.metallicTexPath = main_mat.metallic_texname;
//            main_info.metallic = main_mat.metallic;
//            main_info.roughness = main_mat.roughness;
//            main_info.shininess = main_mat.shininess;
//        }
//        _mat_main = std::make_unique<Material>( dir, "material", main_info );
//
//        //UpdateNormal();
//
//        _vao = std::make_unique<VertexArray>();
//        VertexBufferLayout layout;
//        layout.Push( GL_FLOAT, 3, 0 );
//        layout.Push( GL_FLOAT, 3, sizeof( float ) * 3 );
//        layout.Push( GL_FLOAT, 2, sizeof( float ) * 6 );
//        layout.Push( GL_FLOAT, 3, sizeof( float ) * 8 );
//        layout.Push( GL_FLOAT, 3, sizeof( float ) * 11 );
//        _vbo = std::make_unique<VertexBuffer>( nullptr, 0 );
//        _vao->AddBuffer( *_vbo, layout );
//        UpdateBuffer();
//    }
}

Polyhedron::Polyhedron( const std::string& str )
{
    
}

Polyhedron::Polyhedron( bool render )
    :_render( render )
{
}

Polyhedron::Polyhedron( const Polyhedron& mesh )
    :CPolyhedron( mesh )
{
    _render = mesh._render;
//    if (_render)
//    {
//        MaterialInfos main_info;
//        main_info.diffuseColor = glm::vec3( 255, 205, 196 ) / 255.f;
//        _mat_main = std::make_unique<Material>( "", "material", main_info );
//
//        _vao = std::make_unique<VertexArray>();
//        VertexBufferLayout layout;
//        layout.Push( GL_FLOAT, 3, 0 );
//        layout.Push( GL_FLOAT, 3, sizeof( float ) * 3 );
//        layout.Push( GL_FLOAT, 2, sizeof( float ) * 6 );
//        layout.Push( GL_FLOAT, 3, sizeof( float ) * 8 );
//        layout.Push( GL_FLOAT, 3, sizeof( float ) * 11 );
//        _vbo = std::make_unique<VertexBuffer>( nullptr, 0 );
//        _vao->AddBuffer( *_vbo, layout );
//
//        UpdateNormal();
//        UpdateBuffer();
//    }
}


void Polyhedron::LoadFromData( const float* vertices, const unsigned int* indices, unsigned int nb_vertices, unsigned int nb_faces )
{
    PolyhedronBuilderVF<CPolyhedron::HalfedgeDS> builder{vertices, indices, nb_vertices, nb_faces};
    delegate(builder);
}

void Polyhedron::LoadFromEigen( const Eigen::MatrixXd& V, const Eigen::MatrixXi& F )
{
    PolyhedronBuilderEigen<CPolyhedron::HalfedgeDS> builder(V, F);
    delegate(builder);
}

std::vector<hHalfedge> Polyhedron::FindHoles() const
{
    std::vector<hHalfedge> border_cycles;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(const_cast<Polyhedron&>(*this), std::back_inserter(border_cycles));
    return border_cycles;
}

std::vector<hHalfedge> Polyhedron::HoleEdges( hHalfedge hh ) const
{
    std::vector<hHalfedge> edges;
    if(!hh->is_border())
        return edges;

    for(hHalfedge hh : CGAL::halfedges_around_face(hh, *this))
    {
        edges.push_back(hh);
    }
    return edges;
}

std::pair<std::vector<hFacet>, std::vector<hVertex>> Polyhedron::CloseHole(hHalfedge hh, bool refine, bool refair)
{
    std::vector<hFacet> patch_facets;
    std::vector<hVertex> patch_vertices;
    if(refine && refair)
    {
        CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(*this, hh, std::back_inserter(patch_facets), std::back_inserter(patch_vertices));
    }
    else if (refine)
    {
        CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(*this, hh, std::back_inserter(patch_facets), std::back_inserter(patch_vertices));
    }
    else
    {
        CGAL::Polygon_mesh_processing::triangulate_hole( *this, hh, std::back_inserter( patch_facets ) );
    }
    return std::make_pair(std::move(patch_facets), std::move(patch_vertices));
}

void Polyhedron::LaplacianSmooth( const std::vector<hVertex>& vertices, int nb_iteration, int weight_mode)
{
    for (int i = 0; i < nb_iteration; i++)
    {
        for (auto hv : vertices)
        {
            if (hv->degree() == 0)
                continue;
            float x = 0.f;
            float y = 0.f;
            float z = 0.f;
            float w_sum = 0.f;
            for(auto hh : CGAL::halfedges_around_target(hv->halfedge(), *this))
            {
                auto oppoh = hh->opposite();
                float w = 0.f;
                switch(weight_mode)
                {
                    case 0:
                        w = 1.0f;
                        break;
                    case 1:
                        if(hh->facet() != nullptr)
                            w += FaceArea(hh->facet());
                        if(hh->opposite()->facet() != nullptr)
                            w += FaceArea(hh->opposite()->facet());
                        break;
                    case 2:
                        w = CotWeight(hv, hh);
                        break;
                    default:
                        w = 1.0f;
                }
                x += static_cast<float>(oppoh->vertex()->point().x()) * w;
                y += static_cast<float>(oppoh->vertex()->point().y()) * w;
                z += static_cast<float>(oppoh->vertex()->point().z()) * w;
                w_sum += w;
            }
            x /= w_sum;
            y /= w_sum;
            z /= w_sum;
            hv->point() = { x, y, z };
        }
    }
}

void Polyhedron::LaplacianSmooth( int nb_iteration, int weight_mode )
{
    for (int i = 0; i < nb_iteration; i++)
    {
        for (auto hv = vertices_begin(); hv != vertices_end(); hv++)
        {
            if (hv->degree() == 0)
                continue;

            auto hh = hv->vertex_begin();

            float x = 0.f;
            float y = 0.f;
            float z = 0.f;

            float w_sum = 0.f;
            do
            {
                auto oppoh = hh->opposite();
                if (hh->opposite() == nullptr)
                {
                    continue;
                }
                float w = 0.f;
                switch(weight_mode)
                {
                    case 0:
                        w = 1.0f;
                        break;
                    case 1:
                        w = FaceArea( hh->facet() ) + FaceArea( hh->opposite()->facet() );
                        break;
                    case 2:
                        w = CotWeight(hv, hh);
                        break;
                    default:
                        w = 1.0f;
                }
                x += static_cast<float>(oppoh->vertex()->point().x()) * w;
                y += static_cast<float>(oppoh->vertex()->point().y()) * w;
                z += static_cast<float>(oppoh->vertex()->point().z()) * w;

                w_sum += w;
                hh++;
            } while (hh != hv->vertex_begin());

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

std::array<Polyhedron::Traits::Kernel::Point_3, 8> Polyhedron::ComputeOBB() const
{
    std::array<Traits::Kernel::Point_3, 8> obb_points;
    CGAL::oriented_bounding_box( *this, obb_points, CGAL::parameters::use_convex_hull( true ) );
    return obb_points;
}

CGAL::Bbox_3 Polyhedron::ComputeAABB() const
{
    return CGAL::bbox_3(points_begin(), points_end());
}

std::pair<Polyhedron::Traits::Kernel::Line_3, Polyhedron::Traits::Kernel::Point_3> Polyhedron::ComputePCA() const
{
    Traits::Kernel::Point_3 centroid;
    Traits::Kernel::Line_3 line;
    CGAL::linear_least_squares_fitting_3( points_begin(), points_end(), line, centroid, CGAL::Dimension_tag<0>());
    return {line, centroid};
}

Polyhedron::Traits::Kernel::Point_3 Polyhedron::ComputeCentroidVolumetric() const
{
    using Kernel = Polyhedron::Traits::Kernel;
    Kernel::FT w_sum = 0.0;
    Kernel::Vector_3 pos{ 0.0, 0.0, 0.0 };
    for(auto hf = facets_begin(); hf != facets_end(); hf++)
    {
        Kernel::Point_3 p0 = hf->halfedge()->vertex()->point();
        Kernel::Point_3 p1 = hf->halfedge()->next()->vertex()->point();
        Kernel::Point_3 p2 = hf->halfedge()->prev()->vertex()->point();
        Kernel::Point_3 p3 = CGAL::ORIGIN;

        Kernel::FT w = CGAL::volume(p0, p1, p2, p3);
        w_sum += w;
        pos += w * (p0 - CGAL::ORIGIN);
        pos += w * (p1 - CGAL::ORIGIN);
        pos += w * (p2 - CGAL::ORIGIN);
        pos += w * (p3 - CGAL::ORIGIN);
    }
    pos /= w_sum * 4.0;

    return CGAL::ORIGIN + pos;
}

Polyhedron::Traits::Kernel::Point_3 Polyhedron::ComputeCentroidSurfacial() const
{
    using Kernel = Polyhedron::Traits::Kernel;
    Kernel::FT w_sum = 0.0;
    Kernel::Vector_3 pos{ 0.0, 0.0, 0.0 };
    for(auto hf = facets_begin(); hf != facets_end(); hf++)
    {
        Kernel::Point_3 p0 = hf->halfedge()->vertex()->point();
        Kernel::Point_3 p1 = hf->halfedge()->next()->vertex()->point();
        Kernel::Point_3 p2 = hf->halfedge()->prev()->vertex()->point();

        Kernel::FT w = std::sqrt(CGAL::squared_area(p0, p1, p2));
        w_sum += w;
        pos += w * (p0 - CGAL::ORIGIN);
        pos += w * (p1 - CGAL::ORIGIN);
        pos += w * (p2 - CGAL::ORIGIN);
    }
    pos /= w_sum * 3.0;

    return CGAL::ORIGIN + pos;
}

void Polyhedron::UpdateBuffer()
{
//    if (!_render)
//        return;
//    _vertex_buffer.clear();
//    _vertex_buffer.resize( size_of_facets() * 3 * 14 );
//    int i = 0;
//    for (auto it = facets_begin(); it != facets_end(); it++)
//    {
//        auto hh = it->halfedge();
//        auto hv = hh->vertex();
//        for (int j = 0; j < 3; j++)
//        {
//            _vertex_buffer[i++] = hv->point().x();
//            _vertex_buffer[i++] = hv->point().y();
//            _vertex_buffer[i++] = hv->point().z();
//            _vertex_buffer[i++] = hh->normal.x();
//            _vertex_buffer[i++] = hh->normal.y();
//            _vertex_buffer[i++] = hh->normal.z();
//            _vertex_buffer[i++] = hh->texcoord.x();
//            _vertex_buffer[i++] = hh->texcoord.y();
//            _vertex_buffer[i++] = hv->color.x();
//            _vertex_buffer[i++] = hv->color.y();
//            _vertex_buffer[i++] = hv->color.z();
//            _vertex_buffer[i++] = hh->tangent.x();
//            _vertex_buffer[i++] = hh->tangent.y();
//            _vertex_buffer[i++] = hh->tangent.z();
//
//            hh = hh->next();
//            hv = hh->vertex();
//        }
//    }
//
//    _vbo->UpdateData( _vertex_buffer.data(), _vertex_buffer.size() * sizeof( _vertex_buffer[0] ), GL_STATIC_DRAW );

}


void Polyhedron::UpdateNormal()
{
    CGAL::Polygon_mesh_processing::compute_vertex_normals(*this, VertexNormalMap(this) );
//    for (hVertex hv = vertices_begin(); hv != vertices_end(); hv++)
//    {
//        Halfedge_around_vertex_circulator c( hv->halfedge() );
//        CGAL::Vector_3<Traits::Kernel> normal( 0, 0, 0 );
//        do
//        {
//            if (c == nullptr || c->vertex() == nullptr)
//                continue;
//            Point p0 = c->vertex()->point();
//            Point p1 = c->next()->vertex()->point();
//            Point p2 = c->prev()->vertex()->point();
//            CGAL::Vector_3<Traits::Kernel> p0p1 = p1 - p0;
//            CGAL::Vector_3<Traits::Kernel> p0p2 = p2 - p0;
//            CGAL::Vector_3<Traits::Kernel> face_normal = CGAL::cross_product( p0p1, p0p2 );
//            face_normal /= sqrtf( face_normal.squared_length() );
//            float cos = p0p1.x() * p0p2.x() + p0p1.y() * p0p2.y() + p0p1.z() * p0p2.z();
//            cos = cos / (sqrtf( p0p1.squared_length() * p0p2.squared_length() ));
//            float rad = std::acos( cos );
//            normal += face_normal * rad;
//            c++;
//        } while ((hHalfedge)c != hv->halfedge());
//
//        Eigen::Vector3f n = Eigen::Vector3f( normal.x(), normal.y(), normal.z() );
//        n.normalize();
//        if (n.hasNaN() || n.norm() <= 0.1f)
//        {
//            n = Eigen::Vector3f( 0, 0, 1 );
//        }
//        c = Halfedge_around_vertex_circulator( hv->halfedge() );
//        if (c != nullptr)
//        {
//            do
//            {
//                c->normal = n;
//            } while ((hHalfedge)(++c) != hv->halfedge());
//        }
//    }
    //
    //for (hFacet hf = facets_begin(); hf != facets_end(); hf++)
    //{
    //    hHalfedge hh0 = hf->halfedge();
    //    hHalfedge hh1 = hh0->next();
    //    hHalfedge hh2 = hh1->next();
    //    hVertex hv0 = hh0->vertex();
    //    hVertex hv1 = hh1->vertex();
    //    hVertex hv2 = hh2->vertex();
    //    Point p0 = hv0->point();
    //    Point p1 = hv1->point();
    //    Point p2 = hv2->point();

    //    CGAL::Vector_3<Kernel> face_normal = CGAL::cross_product( p1 - p0, p2 - p0 );

    //    Eigen::Vector3f n( face_normal.x(), face_normal.y(), face_normal.z() );
    //    n.normalize();

    //    hh0->normal = n;
    //    hh1->normal = n;
    //    hh2->normal = n;
    //}

//    for (hFacet hf = facets_begin(); hf != facets_end(); hf++)
//    {
//        hHalfedge hh0 = hf->halfedge();
//        hHalfedge hh1 = hh0->next();
//        hHalfedge hh2 = hh1->next();
//        hVertex hv0 = hh0->vertex();
//        hVertex hv1 = hh1->vertex();
//        hVertex hv2 = hh2->vertex();
//
//        auto temp1 = hv1->point() - hv0->point();
//        auto temp2 = hv2->point() - hv0->point();
//        Eigen::Vector3f edge1( temp1.x(), temp1.y(), temp1.z() );
//        Eigen::Vector3f edge2( temp2.x(), temp2.y(), temp2.z() );
//        Eigen::Vector2f duv1 = hh1->texcoord - hh0->texcoord;
//        Eigen::Vector2f duv2 = hh2->texcoord - hh0->texcoord;
//        float p = 1.0f / (duv1.x() * duv2.y() - duv2.x() * duv1.y());
//        Eigen::Vector3f tangent;
//        tangent.x() = p * (duv2.y() * edge1.x() - duv1.y() * edge2.x());
//        tangent.y() = p * (duv2.y() * edge1.y() - duv1.y() * edge2.y());
//        tangent.z() = p * (duv2.y() * edge1.z() - duv1.y() * edge2.z());
//        tangent.normalize();
//
//        hh0->tangent = tangent;
//        hh1->tangent = tangent;
//        hh2->tangent = tangent;
//    }
}

void Polyhedron::CheckInvalidPoint() const
{
    std::cout << "start checkn\n";
    for (auto it = vertices_begin(); it != vertices_end(); it++)
    {
        if (std::isnan( it->point().x() ) || std::isnan( it->point().y() ) || std::isnan( it->point().z() ))
            std::cout << "HasNan\n";
    }
    std::cout << "check end" << std::endl;
}

void Polyhedron::WriteOFF( const std::string& path ) const
{
#ifndef __EMSCRIPTEN__
    std::ofstream ofs( path );
    CGAL::IO::write_OFF( ofs, *this );
#endif
}

void Polyhedron::WriteOBJ( const std::string& path ) const
{
#ifndef __EMSCRIPTEN__
    std::stringstream ss;

    std::unordered_map<decltype(vertices_begin()), size_t> vertex_id_map;
    size_t v_count = 1;
    for(auto hv = vertices_begin(); hv != vertices_end(); hv++)
    {
        vertex_id_map[hv] = v_count++;
        ss << "v " << hv->point().x() << ' ' << hv->point().y() << ' ' << hv->point().z() << '\n';
    }

    std::unordered_map<decltype(halfedges_begin()), size_t> hedge_id_map;
    size_t h_count = 1;
    for(auto hh = halfedges_begin(); hh != halfedges_end(); hh++)
    {
        hedge_id_map[hh] = h_count++;
        ss << "n " << hh->cgal_n.x() << ' ' << hh->cgal_n.y() << ' ' << hh->cgal_n.z() << '\n';
    }

    for(auto hf = facets_begin(); hf != facets_end(); hf++)
    {
        size_t e0 = hedge_id_map[hf->halfedge()];
        size_t e1 = hedge_id_map[hf->halfedge()->next()];
        size_t e2 = hedge_id_map[hf->halfedge()->prev()];
        size_t v0 = vertex_id_map[hf->halfedge()->vertex()];
        size_t v1 = vertex_id_map[hf->halfedge()->next()->vertex()];
        size_t v2 = vertex_id_map[hf->halfedge()->prev()->vertex()];
        ss << "f " << v0 << "//" << e0 << ' ' << v1 << "//" << e1 << ' ' << v2 << "//" << e2 << '\n';
    }

    std::ofstream ofs(path);
    ofs << ss.rdbuf();
    ofs.close();
#endif
}

int Polyhedron::WriteTriangleSoupVN( float** pvertices, float** pnormals )
{
    int nb_vertices = size_of_facets() * 3;
    *pvertices = new float[nb_vertices * 3];
    *pnormals = new float[nb_vertices * 3];

    int cnt = 0;
    for(auto fh = facets_begin(); fh != facets_end(); fh++)
    {
        Point p[3] = {fh->halfedge()->vertex()->point(), fh->halfedge()->next()->vertex()->point(), fh->halfedge()->prev()->vertex()->point()};
        Traits::Vector_3 n[3] = { fh->halfedge()->cgal_n, fh->halfedge()->next()->cgal_n, fh->halfedge()->prev()->cgal_n };
        for(int i = 0; i < 3; i++)
        {
            (*pvertices)[cnt] = static_cast<float>(p[i].x());
            (*pnormals)[cnt] = n[i].x();
            cnt++;
            (*pvertices)[cnt] = static_cast<float>(p[i].y());
            (*pnormals)[cnt] = n[i].y();
            cnt++;
            (*pvertices)[cnt] = static_cast<float>(p[i].z());
            (*pnormals)[cnt] = n[i].z();
            cnt++;
        }
    }

    return nb_vertices;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> Polyhedron::WriteToEigen() const
{
    Eigen::MatrixXd V(size_of_vertices(), 3);
    int i = 0;
    std::unordered_map<decltype(vertices_begin()), int> idmap;
    for(auto hv = vertices_begin(); hv != vertices_end(); hv++ )
    {
        auto p = hv->point();
        V(i, 0) = p.x();
        V(i, 1) = p.y();
        V(i, 2) = p.z();
        idmap[hv] = i;
        i++;
    }

    Eigen::MatrixXi F(size_of_facets(), 3);
    i = 0;
    for(auto hf = facets_begin(); hf != facets_end(); hf++)
    {
        auto hv0 = hf->halfedge()->vertex();
        auto hv1 = hf->halfedge()->next()->vertex();
        auto hv2 = hf->halfedge()->next()->next()->vertex();
        F(i, 0) = idmap[hv0];
        F(i, 1) = idmap[hv1];
        F(i, 2) = idmap[hv2];
        i++;
    }

    return {std::move(V), std::move(F)};
}

float CotWeight( hVertex vi, hHalfedge hh /*point to vi*/ )
{
    assert( hh->vertex() == vi );

    Eigen::Vector3d pi = ToEigen( hh->vertex()->point() );
    Eigen::Vector3d pj = ToEigen( hh->opposite()->vertex()->point() );
    Eigen::Vector3d pa = ToEigen( hh->next()->vertex()->point() );
    Eigen::Vector3d pb = ToEigen( hh->opposite()->next()->vertex()->point() );

    Eigen::Vector3d papi = pi - pa;
    Eigen::Vector3d papj = pj - pa;

    papi.normalize();
    papj.normalize();

    double cos_a = papi.dot( papj );
    double cot_a = cos_a / (std::sqrt( 1.0 - cos_a * cos_a ));


    Eigen::Vector3d pbpi = pi - pb;
    Eigen::Vector3d pbpj = pj - pb;
    pbpi.normalize();
    pbpj.normalize();
    double cos_b = pbpi.dot( pbpj );
    double cot_b = cos_b / (std::sqrt( 1.0 - cos_b * cos_b ));

    double w = (cot_a + cot_b) / 2.0;
    if (std::isnan( w ) || std::isinf( w ))
        w = 0.f;
    return w;
}

float FaceArea( hFacet hf )
{
    Eigen::Vector3d p0 = ToEigen( hf->halfedge()->vertex()->point() );
    Eigen::Vector3d p1 = ToEigen( hf->halfedge()->next()->vertex()->point() );
    Eigen::Vector3d p2 = ToEigen( hf->halfedge()->prev()->vertex()->point() );
    return std::abs( ((p1 - p0).cross( p2 - p0 )).norm() );
}
