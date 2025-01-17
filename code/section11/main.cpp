
#include "../../external/eigen/Eigen/Dense"
#include <iostream>
 
using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::FullPivLU;


using namespace std;

int main()
{
  // MatrixXd A(3,3) ;
  // MatrixXd b(3,1);
  Matrix3d A;
  Vector3d b;
  A << -12,-10.8,-13.4,
  -18.6,-12.1,-19.6,
  -15.8,-10.4,-11.5;
  b << -0.4,-0.6,-1.4;
  FullPivLU<Matrix3d> lu(A);

  int Arows = A.rows();
  int Acols = A.cols();
  int rankA = lu.rank();

  cout << "Rank of A: " << rankA << endl;

  if (rankA == min(Arows, Acols)){
    Vector3d x = lu.solve(b);
    float error = (A*x-b).cwiseAbs().maxCoeff();
    cout << "Solution Error::" << error << endl;
    cout << "Result: \n" << x << endl;
  }
  else{
    cout  << "Singular matrix A" << endl;
  }
  
  return 0;
}
