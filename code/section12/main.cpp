#include "../../external/eigen/Eigen/Dense"
#include <iostream>

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

int main()
{
  
  Matrix4d C;
  Vector4d d;
  C.row(0) << -2.8,-1.5,-3.6,-1.7;
  C.row(1) << -1.85,-1.2,-1.61,-4.24;
  C.row(2) = 0.4*C.row(0) + 0.23*C.row(1);  
  C.row(3) = -0.5*C.row(0) + 0.2*C.row(2);
  cout << "Matrix C: \n" << C << endl;

  JacobiSVD<Matrix4d, ComputeThinU | ComputeThinV> svd(C);

  d(0) = 5;
  d(1) = -6;
  d(2) = 0.4 * d(0) + 0.23 * d(1);
  d(3) = -0.5 * d(0) + 0.2 * d(2);
  cout << "Vector d: \n" << d << endl;

  FullPivLU<Matrix4d> lu(C);
  int rankC = lu.rank();
  int Crows = C.rows();
  int Ccols = C.cols();

  cout << "Rank of C: " << rankC << endl;

  Vector4d x0 = svd.solve(d);

  float error = (C*x0-d).cwiseAbs().maxCoeff();
  cout << "x0: \n" << x0 << endl;
  cout << "norm of x0: " << x0.norm() << endl;
  cout << "Solution Error:" << error << endl;
  
  MatrixXd nullC = null(C);

  cout << "x0*nullC: \n" << x0.transpose()*nullC << endl;
  return 0;
}

