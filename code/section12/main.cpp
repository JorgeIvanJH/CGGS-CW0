#include "../../external/eigen/Eigen/Dense"
#include <iostream>
#include <cassert>

using namespace Eigen;
using namespace std;

//The mesh quantities
MatrixXi F;
MatrixXd V;

MatrixXd null(const MatrixXd& A){
  JacobiSVD<Matrix4d, ComputeThinU | ComputeThinV> svd(A); // SVD Decomposition
  double threshold = 1e-10;
  Eigen::VectorXd singularValues = svd.singularValues();
  Eigen::MatrixXd V = svd.matrixV();
  Eigen::MatrixXd nullA = V.block(0, svd.rank(), V.rows(), V.cols() - svd.rank());
  return nullA;
  
}

double error_calc(const MatrixXd& A){
  return A.cwiseAbs().maxCoeff();
}

int main()
{
  
  Matrix4d C;
  Vector4d d;

  // Left hand side matrix
  C.row(0) << -2.8,-1.5,-3.6,-1.7;
  C.row(1) << -1.85,-1.2,-1.61,-4.24;
  C.row(2) = 0.4*C.row(0) + 0.23*C.row(1);  
  C.row(3) = -0.5*C.row(0) + 0.2*C.row(2);
  cout << "Matrix C: " << C << endl;
  FullPivLU<Matrix4d> lu(C);
  int rankC = lu.rank();
  int Crows = C.rows();
  int Ccols = C.cols();
  cout << "Rank of C: " << rankC << endl;

  // Right hand side vector
  d(0) = 5;
  d(1) = -6;
  d(2) = 0.4 * d(0) + 0.23 * d(1);
  d(3) = -0.5 * d(0) + 0.2 * d(2);
  cout << "Vector d: " << d << endl;

  // First Solution Using LU Decomposition
  cout << "LU Decomposition" << endl;
  Vector4d x_lu = lu.solve(d);
  cout << "x_lu: " << x_lu << endl;
  cout << "norm of x_lu: " << x_lu.norm() << endl;

  // Second Solution Using SVD Decomposition
  cout << "SVD Decomposition" << endl;
  JacobiSVD<Matrix4d, ComputeThinU | ComputeThinV> svd(C);
  Vector4d x_svd = svd.solve(d);
  cout << "x_svd: " << x_svd << endl;
  cout << "norm of x_svd: " << x_svd.norm() << endl;
  float error = error_calc(C*x_svd-d);
  cout << "Solution Error: " << error << endl;
  MatrixXd nullC = null(C);
  cout << "x0*nullC: " << x_svd.transpose()*nullC << endl;

  // Validating with the null space
  std::array<double, 3> alphas = {0.5, -0.3, -3.981};
  std::array<double, 3> betas = {0.5, 10.2, 6.1};
  for (int i = 0; i < 3; i++){
    cout << "(alpha, beta): (" << alphas[i] << ", " << betas[i] << ")" << endl;
    Vector2d alpha_beta;
    alpha_beta << alphas[i], betas[i];
    Vector4d x = x_svd + nullC*alpha_beta;
    error = error_calc(C*x-d);
    cout << "error: " << error << endl;
    assert(error < 10e-10 && "error is not neglible");
  }

  // Incompatible Right hand side vector
  cout << "Incompatible Right hand side vector" << endl;
  d(3) = 0.3945;
  x_svd = svd.solve(d);
  cout << "x_svd: " << x_svd << endl;
  error = error_calc(C*x_svd-d);
  cout << "error: " << error << endl;
  double enorm ; 
  enorm = (C.transpose()*C*x_svd - C.transpose()*d).norm();
  cout << "C.transpose()*C*xe - C.transpose()*d error: " << enorm << endl;
  return 0;
}

