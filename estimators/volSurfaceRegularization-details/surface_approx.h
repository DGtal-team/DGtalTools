#pragma once

#include <DGtal/helpers/StdDefs.h>
#include <DGtal/topology/CanonicCellEmbedder.h>
#include <DGtal/topology/CanonicSCellEmbedder.h>
#include <DGtal/math/linalg/EigenSupport.h>
#include <DGtal/dec/DiscreteExteriorCalculus.h>
#include <tuple>

bool ends_with(const std::string& value, const std::string& ending);

typedef DGtal::DiscreteExteriorCalculus<2, 3, DGtal::EigenLinearAlgebraBackend> Calculus;
typedef Calculus::LinearAlgebraBackend Backend;
typedef Calculus::KSpace KSpace;
typedef DGtal::CanonicCellEmbedder<KSpace> CellEmbedder;
typedef DGtal::CanonicSCellEmbedder<KSpace> SCellEmbedder;
typedef KSpace::SCell SCell;
typedef KSpace::Cell Cell;
typedef Backend::DenseVector FlatVector;
typedef Backend::SparseMatrix OperatorMatrix;
typedef Backend::Triplet Triplet;
typedef DGtal::Z3i::Point Point;
typedef DGtal::Z3i::RealPoint RealPoint;

std::tuple<Calculus, FlatVector>
initCalculusAndNormalsFromSurfelNormalsCSV(const std::string& filename);

bool
checkOperatorSymmetry(const OperatorMatrix& matrix, const double tol=1e-8);

FlatVector
vertexNormals(const Calculus& calculus, const FlatVector& face_normals);

struct ApproxParams
{
    double regularization_position;
    double regularization_center;
    double align;
    double fairness;
    double barycenter;
};

void
exportOBJ(const Calculus& calculus, const FlatVector& positions, const std::string& filename);

// (original_positions, regularization_position, original_centers, regularized_centers)
std::tuple<FlatVector, FlatVector, FlatVector, FlatVector>
approximateSurface(const Calculus& calculus, const FlatVector& normals, const ApproxParams& params);

// (L_position, L_align)
std::tuple<double, double>
approximateSurfaceEnergies(const Calculus& calculus, const FlatVector& normals, const FlatVector& positions);

