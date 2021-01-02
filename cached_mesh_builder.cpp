/**
 * @file    cached_mesh_builder.cpp
 *
 * @author  Martin Hošala <xhosal00@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using pre-computed field
 *
 * @date    15. 12. 2020
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "cached_mesh_builder.h"

CachedMeshBuilder::CachedMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Cached")
{

}

void CachedMeshBuilder::precalculateFieldDistances(const ParametricScalarField &field)
{
  const Vec3_t<float> *pPoints = field.getPoints().data();
  const unsigned count = unsigned(field.getPoints().size());
  mFieldDistances.resize((mGridSize+1)*(mGridSize+1)*(mGridSize+1));
  #pragma omp parallel for
  for (size_t x = 0; x <= mGridSize; x++) {
    for (size_t y = 0; y <= mGridSize; y++) {
      for (size_t z = 0; z <= mGridSize; z++) {
        float value = std::numeric_limits<float>::max();
        for(unsigned i = 0; i < count; ++i)
        {
          float distanceSquared  = (x * mGridResolution - pPoints[i].x) * (x * mGridResolution - pPoints[i].x);
          distanceSquared       += (y * mGridResolution - pPoints[i].y) * (y * mGridResolution - pPoints[i].y);
          distanceSquared       += (z * mGridResolution - pPoints[i].z) * (z * mGridResolution - pPoints[i].z);

          value = std::min(value, distanceSquared);
        }
        mFieldDistances[x*(mGridSize+1)*(mGridSize+1) + y*(mGridSize+1) + z] = sqrt(value);
      }
    }
  }
}

unsigned CachedMeshBuilder::marchCubes(const ParametricScalarField &field)
{
  precalculateFieldDistances(field);

  size_t totalCubesCount = mGridSize*mGridSize*mGridSize;
  unsigned totalTriangles = 0;

  #pragma omp parallel for reduction(+:totalTriangles)
  for(size_t i = 0; i < totalCubesCount; ++i)
  {
      Vec3_t<float> cubeOffset( i % mGridSize,
                               (i / mGridSize) % mGridSize,
                                i / (mGridSize*mGridSize));

      totalTriangles += buildCube(cubeOffset, field);
  }

  return totalTriangles;
}

float CachedMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
  unsigned x = std::floor(pos.x/mGridResolution + 1/2);
  unsigned y = std::floor(pos.y/mGridResolution + 1/2);
  unsigned z = std::floor(pos.z/mGridResolution + 1/2);
  return mFieldDistances[x*(mGridSize+1)*(mGridSize+1) + y*(mGridSize+1) + z];
}

void CachedMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
  #pragma omp critical
  {
    mTriangles.push_back(triangle);
  }
}
