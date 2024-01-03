//
// Created by yrz on 3/19/23.
//

#ifndef ELASTICITY_GUMGEN_SUPPORT_CLASSES_H
#define ELASTICITY_GUMGEN_SUPPORT_CLASSES_H
#include "Polyhedron.h"
#include <tuple>
#include <Eigen/Eigen>

class OralScanMesh : public Polyhedron {
public:
    friend class FakeGumGen;
    OralScanMesh( bool upper );
    void ReadJsonLabels( const std::string& path );
    void ReadLabels(const unsigned int* labels);
    void VertexLabelToFaceLabel();
    void ComputeOrientation();
    void AlignMesh();
    void AlignMeshInverse();
    //Remove gum and wisdom teeth
    void RemoveGumPart();

private:
    Eigen::Vector3d _obb_center;
    Eigen::Vector3d _obb_up;
    bool _upper;
};

#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME OralScanMesh
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME Polyhedron
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>
#undef CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME
#undef CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME

class TeethMeshes
{
public:
    TeethMeshes(bool upper);
    void InitFromOralScanMesh(OralScanMesh& oral_scan_mesh);
    void InitFromFile( const std::array<std::string, 49>& paths);
    void InitFromString( const std::array<std::string, 49>& strings);
    void AlignMeshes( Eigen::Vector3d obb_up, Eigen::Vector3d obb_center);
    void AlignMeshesInverse();
    void ComputeTeethOBB();
    void LoadTeethOritFromJson( const char* str );
    void ComputeTeethPCA();
    std::tuple<Eigen::Vector3f, Eigen::Vector3f, Eigen::Vector3f> ComputeUpdirAndCenterPointsFromOBB( const std::array<Eigen::Vector3f, 8>& obb);

    std::array<std::array<Eigen::Vector3f, 8>, 49> _teeth_obb;
    std::array<Eigen::Vector3f, 49> _teeth_updir;
    std::array<Eigen::Vector3f, 49> _teeth_leftdir;
    std::array<Eigen::Vector3f, 49> _teeth_bottom_centers;
    std::array<Eigen::Vector3f, 49> _teeth_pca_dir;
    std::array<Eigen::Vector3f, 49> _teeth_pca_centroid;
    std::array<std::unique_ptr<Polyhedron>, 49> _teeth;
    Eigen::Vector3d _obb_up;
    Eigen::Vector3d _obb_center;
    bool _upper;
};
#endif //ELASTICITY_GUMGEN_SUPPORT_CLASSES_H