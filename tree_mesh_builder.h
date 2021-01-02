/**
 * @file    tree_mesh_builder.h
 *
 * @author  Martin Hošala <xhosal00@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    14. 12. 2020
 **/

#ifndef TREE_MESH_BUILDER_H
#define TREE_MESH_BUILDER_H

#include <vector>
#include "base_mesh_builder.h"

class TreeMeshBuilder : public BaseMeshBuilder
{
public:
    TreeMeshBuilder(unsigned gridEdgeSize);

protected:
    unsigned marchCubes(const ParametricScalarField &field);
    float evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field);
    void emitTriangle(const Triangle_t &triangle);
    const Triangle_t *getTrianglesArray() const { return mTriangles.data(); }

    std::vector<Triangle_t> mTriangles; ///< Temporary array of triangles
    unsigned handleNode(const ParametricScalarField &field, const Vec3_t<float> &nodePosition, const unsigned &nodeEdgeSize);
    unsigned mCutOff;
};

#endif // TREE_MESH_BUILDER_H
