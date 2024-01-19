//
// Created by yrz on 7/7/22.
//

#ifndef ELASTICITY_POLYHEDRON_H
#define ELASTICITY_POLYHEDRON_H
#include <memory>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <Eigen/Eigen>
#include <Eigen/Core>

using Real = double;

template <typename SizeType>
struct TTriangle
{
public:
    using size_type = SizeType;
    TTriangle(SizeType i0, SizeType i1, SizeType i2)
    {
        _id[0] = i0;
        _id[1] = i1;
        _id[2] = i2;
    }
    SizeType &operator[](uint32_t i) { return _id[i]; }
    const SizeType &operator[](uint32_t i) const { return _id[i]; }
    std::pair<SizeType, SizeType> GetEdge(uint32_t i)
    {
        switch (i)
        {
        case 0:
            return std::make_pair(_id[0], _id[1]);
            break;
        case 1:
            return std::make_pair(_id[1], _id[2]);
            break;
        case 2:
            return std::make_pair(_id[2], _id[0]);
            break;
        }
        return std::make_pair<size_t, size_t>(0, 0);
    }

protected:
    SizeType _id[3]{0, 0, 0};
};
template <typename SizeType>
struct TEdge
{
public:
    using size_type = SizeType;

    TEdge() = default;

    TEdge(SizeType i0, SizeType i1)
    {
        _i0 = i0;
        _i1 = i1;
    }

    SizeType _i0;
    SizeType _i1;
    std::vector<SizeType> _faces;
};

template <typename SizeType>
struct TPairHashUnordered
{
    using size_type = SizeType;
    size_t operator()(const std::pair<SizeType, SizeType> &p) const
    {
        if (p.first > p.second)
        {
            return boost::hash<std::pair<SizeType, SizeType>>()(p);
        }
        else
        {
            return boost::hash<std::pair<SizeType, SizeType>>()({p.second, p.first});
        }
    }
};

template <typename SizeType>
struct TPairPredUnordered
{
    using size_type = SizeType;

    bool operator()(std::pair<SizeType, SizeType> e0, std::pair<SizeType, SizeType> e1) const
    {
        return e0.first == e1.first && e0.second == e1.second || e0.first == e1.second && e0.second == e1.first;
    }
};

template <typename SizeType>
struct TPairHash
{
    using size_type = SizeType;
    size_t operator()(std::pair<SizeType, SizeType> edge) const
    {
        return boost::hash<std::pair<SizeType, SizeType>>()(edge);
    }
};

template <typename SizeType>
struct TPairPred
{
    using size_type = SizeType;
    bool operator()(std::pair<SizeType, SizeType> edge0, std::pair<SizeType, SizeType> edge1) const
    {
        return edge0.first == edge1.first && edge0.second == edge1.second;
    }
};

template <class Refs, typename Tag, typename Point>
class MyVertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>
{
public:
    MyVertex() = default;
    explicit MyVertex( const Point& p ) : CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>( p ) {}

public:
    int label{ 0 };
};

template<class Refs, class Traits>
class MyHalfedge : public CGAL::HalfedgeDS_halfedge_base<Refs>
{
public:
    typename Traits::Vector_3 cgal_n;
};

template<class Refs, class Traits>
class MyFace : public CGAL::HalfedgeDS_face_base<Refs>
{
public:
    int label{ 0 };
};

class MyItems : public CGAL::Polyhedron_items_3
{
public:
    template<class Refs, class Traits>
    struct Vertex_wrapper
    {
        using Point = typename Traits::Point_3;
        using Vertex = MyVertex<Refs, CGAL::Tag_true, Point>;
    };
    template <typename Refs, typename Traits>
    struct Halfedge_wrapper
    {
        using Halfedge = MyHalfedge<Refs, Traits>;
    };
    template <typename Refs, typename Traits>
    struct Face_wrapper
    {
        using Face = MyFace<Refs, Traits>;
    };

};
using KernelSimpleCartesian = CGAL::Simple_cartesian<Real>;
using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Kernel = KernelEpick;
using CPolyhedron = CGAL::Polyhedron_3<KernelEpick, MyItems>;
using hHalfedge = CPolyhedron::Halfedge_handle;
using hVertex = CPolyhedron::Vertex_handle;
using hCVertex = CPolyhedron::Vertex_const_handle;
using hFacet = CPolyhedron::Facet_handle;
using Halfedge = CPolyhedron::Halfedge;
using CVertex = CPolyhedron::Vertex;
using Facet = CPolyhedron::Facet;
using iHalfedge = CPolyhedron::Halfedge_iterator;
using iVertex = CPolyhedron::Vertex_iterator;
using iFacet = CPolyhedron::Facet_iterator;
using Point3 = CPolyhedron::Point_3;
using Vec3 = CPolyhedron::Traits::Vector_3;

struct VertexNormalMap
{
    VertexNormalMap(CPolyhedron* mesh) : _mesh(mesh)
    {}
    using value_type = CPolyhedron::Traits::Vector_3;
    using key_type = boost::graph_traits<CPolyhedron>::vertex_descriptor;
    CPolyhedron* _mesh;
};

inline void put(VertexNormalMap m, VertexNormalMap::key_type hv, VertexNormalMap::value_type n)
{
    for(auto hh : CGAL::halfedges_around_target(hv, *(m._mesh)))
    {
        hh->cgal_n = n;
    }
}
template <typename HDS, typename Kernel>
class TPolyhedronObjBuilder : public CGAL::Modifier_base<HDS>
{
public:
    TPolyhedronObjBuilder(const std::vector<typename Kernel::Point_3> &vertices, const std::vector<size_t> &indices)
        : _vertices(vertices), _indices(indices) {}
    virtual void operator()(HDS &hds) override
    {
        _success = true;
        CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
        builder.begin_surface(_vertices.size(), _indices.size() / 3);
        for (size_t i = 0, size = _vertices.size(); i < size; i += 1)
        {
            if(builder.add_vertex(_vertices[i]) == nullptr)
            {
                _success = false;
                printf("Error: failed to add vertex: %zd\n", i);
                break;
            }
        }
        if(!_success)
        {
            return;
        }
        for (size_t f = 0, size = _indices.size() / 3; f < size; ++f)
        {
            builder.begin_facet();
            builder.add_vertex_to_facet(_indices[f * 3 + 0]);
            builder.add_vertex_to_facet(_indices[f * 3 + 1]);
            builder.add_vertex_to_facet(_indices[f * 3 + 2]);
            typename HDS::Halfedge_handle hh = builder.end_facet();
            if (hh == nullptr)
            {
                _success = false;
                printf("Error: failed to add face: %zd\n", f);
                break;
            }
        }
        builder.end_surface();
        if(builder.error())
        {
            _success = false;
        }
    }

    bool Success() const { return _success; };

protected:
    const std::vector<typename Kernel::Point_3> &_vertices;
    const std::vector<size_t> &_indices;
    bool _success = false;
};

class Polyhedron : public CPolyhedron
{
public:
    using Triangle = TTriangle<typename CPolyhedron::Vertex::size_type>;
    using PolyhedronObjBuilder = TPolyhedronObjBuilder<typename CPolyhedron::HalfedgeDS, Traits::Kernel>;
    Polyhedron( const std::string& path, bool render );
    explicit Polyhedron( const std::string& str );
    explicit Polyhedron( bool render = false);
    Polyhedron( const Polyhedron& mesh );

    void LoadFromData( const float* vertices, const unsigned int* indices, unsigned int nb_vertices, unsigned int nb_faces );
    void LoadFromEigen( const Eigen::MatrixXd& V, const Eigen::MatrixXi& F );
    void BuildFromVerticesIndices(const std::vector<Traits::Kernel::Point_3> &vertices, const std::vector<typename CPolyhedron::Vertex::size_type> &indices)
    {
        this->clear();
        PolyhedronObjBuilder builder(vertices, indices);
        this->delegate(builder);
        if(!builder.Success())
        {
            
        }
    }
    void BuildFromVerticesFaces(const std::vector<Traits::Kernel::Point_3>& vertices, const std::vector<Triangle>& faces)
    {
        this->clear();
        std::vector<typename CPolyhedron::Vertex::size_type> indices;
        for(const auto& f : faces)
        {
            indices.push_back(f[0]);
            indices.push_back(f[1]);
            indices.push_back(f[2]);
        }
        BuildFromVerticesIndices(vertices, indices);
    }
    void LoadLabels(const std::vector<int>& labels);
    std::vector<hHalfedge> FindHoles() const;
    std::vector<hHalfedge> HoleEdges(hHalfedge hh ) const;
    std::pair<std::vector<hFacet>, std::vector<hVertex>> CloseHole(hHalfedge hh, bool refine = false, bool refair = false);
    void LaplacianSmooth( const std::vector<hVertex>& vertices, int nb_iteration, int weight_mode);
    void LaplacianSmooth( int nb_iteration, int weight_mode );

    std::array<Point_3, 8> ComputeOBB() const;
    CGAL::Bbox_3 ComputeAABB() const;
    std::pair<Traits::Kernel::Line_3, Traits::Kernel::Point_3> ComputePCA() const;
    Traits::Kernel::Point_3 ComputeCentroidVolumetric() const;
    Traits::Kernel::Point_3 ComputeCentroidSurfacial() const;
    bool IsSmallHole(typename CPolyhedron::Halfedge_handle hh, int max_num_hole_edges, float max_hole_diam)
    {
        int num_hole_edges = 0;
        CGAL::Bbox_3 hole_bbox;
        for (typename CPolyhedron::Halfedge_handle hc : CGAL::halfedges_around_face(hh, (CPolyhedron&)*this))
        {
            const typename Kernel::Point_3& p = hc->vertex()->point();
            hole_bbox += p.bbox();
            ++num_hole_edges;
            // Exit early, to avoid unnecessary traversal of large holes
            if (num_hole_edges > max_num_hole_edges) return false;
            if (hole_bbox.xmax() - hole_bbox.xmin() > max_hole_diam) return false;
            if (hole_bbox.ymax() - hole_bbox.ymin() > max_hole_diam) return false;
            if (hole_bbox.zmax() - hole_bbox.zmin() > max_hole_diam) return false;
        }
        return true;
    }
    void UpdateBuffer();
    void UpdateNormal();
    void CheckInvalidPoint() const;
    void WriteOFF( const std::string& path ) const;
    void WriteOBJ( const std::string& path ) const;
    int WriteTriangleSoupVN( float** vertices, float** normals );
    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> WriteToEigen() const;
    std::pair<std::vector<Traits::Kernel::Point_3>, std::vector<Triangle>> ToVerticesTriangles() const
    {
        std::unordered_map<hCVertex, CPolyhedron::Vertex::size_type> idmap;
        std::vector<Traits::Kernel::Point_3> vertices;
        std::vector<Triangle> triangles;
        size_t cnt = 0;
        for (auto hv = this->vertices_begin(); hv != this->vertices_end(); hv++)
        {
            vertices.push_back(hv->point());
            idmap.insert({(hCVertex)hv, cnt++});
        }
        for (auto hf = this->facets_begin(); hf != this->facets_end(); hf++)
        {
            triangles.emplace_back(
                idmap.at((hCVertex)(hf->halfedge()->vertex())),
                idmap.at((hCVertex)(hf->halfedge()->next()->vertex())),
                idmap.at((hCVertex)(hf->halfedge()->prev()->vertex())));
        }
        return {vertices, triangles};
    }


private:
    bool _render;
//    std::unique_ptr<VertexArray> _vao{ nullptr };
//    std::unique_ptr<VertexBuffer> _vbo{ nullptr };
    std::vector<float> _vertex_buffer;
};

float CotWeight( hVertex vi, hHalfedge hh );

float FaceArea( hFacet hf );


template <typename HDS>
class PolyhedronBuilderVF : public CGAL::Modifier_base<HDS>
{
public:
    PolyhedronBuilderVF( const float* vertices, const unsigned int* indices, unsigned int nb_vertices, unsigned int nb_faces)
        :_vertices(vertices), _indices(indices), _nb_vertices(nb_vertices), _nb_faces(nb_faces) {}
    virtual void operator()(HDS& hds) override;
private:
    const float* _vertices{nullptr};
    const unsigned int* _indices{nullptr};
    unsigned int _nb_vertices{0};
    unsigned int _nb_faces{0};
};

template <typename HDS>
inline void PolyhedronBuilderVF<HDS>::operator()(HDS& hds)
{
    CGAL::Polyhedron_incremental_builder_3<HDS> builder( hds, true );
    builder.begin_surface( _nb_vertices, _nb_faces );

    for (size_t i = 0; i < _nb_vertices; i++)
    {
        hVertex vh = builder.add_vertex( Point3( _vertices[i * 3 + 0], _vertices[i * 3 + 1], _vertices[i * 3 + 2] ) );
    }

    for (unsigned f = 0; f < _nb_faces; ++f)
    {
        builder.begin_facet();
        builder.add_vertex_to_facet( _indices[f * 3 + 0] );
        builder.add_vertex_to_facet( _indices[f * 3 + 1] );
        builder.add_vertex_to_facet( _indices[f * 3 + 2] );
        hHalfedge hh = builder.end_facet();
    }
    builder.end_surface();
}

template <typename HDS>
class PolyhedronBuilderEigen : public CGAL::Modifier_base<HDS>
{
public:
    PolyhedronBuilderEigen( const Eigen::MatrixXd& V, const Eigen::MatrixXi& F );
    virtual void operator()(HDS& hds) override;
protected:
    const Eigen::MatrixXd& _V;
    const Eigen::MatrixXi& _F;
};

template <typename HDS>
inline PolyhedronBuilderEigen<HDS>::PolyhedronBuilderEigen( const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) 
    :_V(V), _F(F) {}

template <typename HDS>
inline void PolyhedronBuilderEigen<HDS>::operator()(HDS& hds)
{
    CGAL::Polyhedron_incremental_builder_3<HDS> builder( hds, true );
    builder.begin_surface(_V.rows(), _F.rows());

    for(Eigen::Index i = 0; i < _V.rows(); i++)
    {
        auto vh = builder.add_vertex(Point3( _V(i, 0), _V(i, 1), _V(i, 2) ));
    }

    for(Eigen::Index i = 0; i < _F.rows(); i++)
    {
        builder.begin_facet();
        builder.add_vertex_to_facet( _F(i, 0) );
        builder.add_vertex_to_facet( _F(i, 1) );
        builder.add_vertex_to_facet( _F(i, 2) );
        builder.end_facet();
    }
    builder.end_surface();
}


#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME Polyhedron
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CPolyhedron
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>
#undef CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME
#undef CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME
#endif //ELASTICITY_POLYHEDRON_H
