#ifndef HARMONIC_INTERPOLATION_HEADER_FILE
#define HARMONIC_INTERPOLATION_HEADER_FILE

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "set_diff.h"
#include "slice_columns_sparse.h"

// Helper function to construct the differential matrix d0
Eigen::SparseMatrix<double> construct_d0(const Eigen::MatrixXi &E, int numVertices)
{
  std::vector<Eigen::Triplet<double>> triplets; // triples to construct the sparse matrix
  for (int i = 0; i < E.rows(); ++i)
  {
    int source = E(i, 0);                   // Source vertex
    int target = E(i, 1);                   // Target vertex
    triplets.emplace_back(i, source, -1.0); // d0(i, source) = -1
    triplets.emplace_back(i, target, 1.0);  // d0(i, target) = 1
  }
  Eigen::SparseMatrix<double> d0(E.rows(), numVertices);
  d0.setFromTriplets(triplets.begin(), triplets.end());
  return d0;
}

Eigen::VectorXd harmonic_interpolation(const Eigen::MatrixXd &V,
                                       const Eigen::MatrixXi &E,
                                       const Eigen::VectorXi &B,
                                       const Eigen::VectorXd &xB)
{
  int numVertices = V.rows();

  // Step 1: Compute I (complement of B)
  Eigen::VectorXi I = set_diff(Eigen::VectorXi::LinSpaced(numVertices, 0, numVertices - 1), B);

  // Step 2: Construct the differential matrix d0
  Eigen::SparseMatrix<double> d0 = construct_d0(E, numVertices);

  // Step 3: Slice d0 into d0|I and d0|B
  Eigen::SparseMatrix<double> d0I = slice_columns_sparse(d0, I);
  Eigen::SparseMatrix<double> d0B = slice_columns_sparse(d0, B);

  // Step 4: Compute the RHS of the equation
  Eigen::VectorXd rhs = -d0I.transpose() * (d0B * xB);

  // Step 5: Solve the system (d0I^T * d0I) * xI = rhs
  Eigen::SparseMatrix<double> A = d0I.transpose() * d0I;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(A);
  Eigen::VectorXd xI = solver.solve(rhs);

  // Step 6: Combine xB and xI into the final result x
  Eigen::VectorXd x(numVertices);
  for (int i = 0; i < B.size(); ++i)
  {
    x[B[i]] = xB[i]; // Set fixed values
  }
  for (int i = 0; i < I.size(); ++i)
  {
    x[I[i]] = xI[i]; // Set interpolated values
  }

  return x;
}

#endif
