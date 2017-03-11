#pragma once

#include <tuple>

template <typename Shape>
Calculus
createCalculusFromShapeBorder(const KSpace& kspace, const Shape& shape);

template <typename Shape>
FlatVector
computeFaceNormals(const Calculus& calculus, const Shape& shape, const double radius);

template <typename Shape>
std::tuple<Calculus, FlatVector>
initCalculusAndNormals(const KSpace& kspace, const Shape& shape, const double radius);

template <typename Shape>
std::tuple<Calculus, FlatVector>
initCalculusAndNormalsWithNoise(const KSpace& kspace, const Shape& shape, const double radius, const double noise);

#include "surface_extract.ih"

