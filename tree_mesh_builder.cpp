/**
 * @file    tree_mesh_builder.cpp
 *
 * @author  Martin Hošala <xhosal00@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    14. 12. 2020
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "tree_mesh_builder.h"

TreeMeshBuilder::TreeMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Octree")
{
  mCutOff = 1;
}
unsigned TreeMeshBuilder::handleNode(const ParametricScalarField &field, const Vec3_t<float> &nodePosition, const unsigned &nodeEdgeSize)
{
  if (nodeEdgeSize <= mCutOff)
    return buildCube(nodePosition, field);

  auto halfEdge = nodeEdgeSize / 2;
  auto x = nodePosition.x, y = nodePosition.y, z = nodePosition.z;
  auto xCenter = x + halfEdge, yCenter = y + halfEdge, zCenter = z + halfEdge;
  auto nodeCenter = Vec3_t<float>(xCenter * mGridResolution,
                                  yCenter * mGridResolution,
                                  zCenter * mGridResolution);

  auto F_p = evaluateFieldAt(nodeCenter, field);

  if (F_p > mIsoLevel + sqrt(3)/2 * nodeEdgeSize * mGridResolution)
    return 0;

  Vec3_t<float> childNodes[] = {
    Vec3_t<float>(x      , y      , z      ),
    Vec3_t<float>(x      , y      , zCenter),
    Vec3_t<float>(x      , yCenter, z      ),
    Vec3_t<float>(x      , yCenter, zCenter),
    Vec3_t<float>(xCenter, y      , z      ),
    Vec3_t<float>(xCenter, y      , zCenter),
    Vec3_t<float>(xCenter, yCenter, z      ),
    Vec3_t<float>(xCenter, yCenter, zCenter)
  };

  unsigned childrenSums[8] = {0};

  for (size_t i = 0; i < 8; i++) {
    #pragma omp task shared(childrenSums)
    childrenSums[i] += handleNode(field, childNodes[i], halfEdge);
  }
  #pragma omp taskwait

  unsigned childrenSum = 0;
  for (size_t i = 0; i < 8; i++)
      childrenSum += childrenSums[i];

  return childrenSum;
}

unsigned TreeMeshBuilder::marchCubes(const ParametricScalarField &field)
{
    // Suggested approach to tackle this problem is to add new method to
    // this class. This method will call itself to process the children.
    // It is also strongly suggested to first implement Octree as sequential
    // code and only when that works add OpenMP tasks to achieve parallelism.
    // 1. Compute total number of cubes in the grid.
    // size_t totalCubesCount = mGridSize*mGridSize*mGridSize;
    //
    // unsigned totalTriangles = 0;
    //
    // // 2. Loop over each coordinate in the 3D grid.
    // for(size_t i = 0; i < totalCubesCount; ++i)
    // {
    //     // 3. Compute 3D position in the grid.
    //     Vec3_t<float> cubeOffset( i % mGridSize,
    //                              (i / mGridSize) % mGridSize,
    //                               i / (mGridSize*mGridSize));
    //
    //     // 4. Evaluate "Marching Cube" at given position in the grid and
    //     //    store the number of triangles generated.
    //     totalTriangles += buildCube(cubeOffset, field);
    // }

    // 5. Return total number of triangles generated.
    unsigned count = 0;
    # pragma omp parallel
    {
        # pragma omp single
        {
          count = handleNode(field, {0, 0, 0}, mGridSize);
        }
    }
    return count;
}

float TreeMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
  // NOTE: This method is called from "buildCube(...)"!

  // 1. Store pointer to and number of 3D points in the field
  //    (to avoid "data()" and "size()" call in the loop).
  const Vec3_t<float> *pPoints = field.getPoints().data();
  const unsigned count = unsigned(field.getPoints().size());

  float value = std::numeric_limits<float>::max();

  // 2. Find minimum square distance from points "pos" to any point in the
  //    field.
  for(unsigned i = 0; i < count; ++i)
  {
      float distanceSquared  = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
      distanceSquared       += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
      distanceSquared       += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

      // Comparing squares instead of real distance to avoid unnecessary
      // "sqrt"s in the loop.
      value = std::min(value, distanceSquared);
  }

  // 3. Finally take square root of the minimal square distance to get the real distance
  return sqrt(value);
}

void TreeMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
  #pragma omp critical
  {
    mTriangles.push_back(triangle);
  }
}
