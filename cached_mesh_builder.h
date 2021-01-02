/**
 * @file    cached_mesh_builder.h
 *
 * @author  Martin Hošala <xhosal00@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using pre-computed field
 *
 * @date    15. 12. 2020
 **/

#ifndef CACHED_MESH_BUILDER_H
#define CACHED_MESH_BUILDER_H

#include <vector>
#include "base_mesh_builder.h"

class CachedMeshBuilder : public BaseMeshBuilder
{
public:
    CachedMeshBuilder(unsigned gridEdgeSize);

protected:
    unsigned marchCubes(const ParametricScalarField &field);
    float evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field);
    void emitTriangle(const Triangle_t &triangle);
    const Triangle_t *getTrianglesArray() const { return mTriangles.data(); }

    std::vector<Triangle_t> mTriangles; ///< Temporary array of triangles
    std::vector<float> mFieldDistances; ///< Temporary array of triangles

    void precalculateFieldDistances(const ParametricScalarField &field);
};

#endif // CACHED_MESH_BUILDER_H
