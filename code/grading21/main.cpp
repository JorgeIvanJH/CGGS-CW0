#include "readOFF.h"
#include "compute_areas_normals.h"
#include "serialization.h"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <filesystem>

using namespace Eigen;
using namespace std;

//The mesh quantities
MatrixXi F;
MatrixXd V;
MatrixXd normals;
VectorXd areas;

double tolerance = 1e-8;

namespace fs = std::filesystem;

int main()
{
  
  std::string folderPath(DATA_PATH); // Replace with your folder path
  for (const auto& entry : fs::directory_iterator(folderPath)) {
    if (entry.is_regular_file() && entry.path().extension() == ".off") {
      cout<<"Working on file "<<entry.path().filename()<<endl;
      std::string dataName = entry.path().string();
      dataName.erase(dataName.size() - 4, 4);
      std::ifstream ifs(dataName+"-section21.data", std::ofstream::binary);
      MatrixXd normalsGT;
      VectorXd areasGT;
      readOFF(entry.path().string(), V, F);
      compute_areas_normals(V, F, areas, normals);
      deserializeVector(areasGT, ifs);
      deserializeMatrix(normalsGT, ifs);
      ifs.close();
      
      //comparisons and feedback
      if ((areas.rows()!=areasGT.rows())||(areas.cols()!=areasGT.cols())){
        cout<<"areas is of the wrong size!"<<endl;
        continue;
      }
      if ((normals.rows()!=normalsGT.rows())||(normals.cols()!=normalsGT.cols())){
        cout<<"normals is of the wrong size!"<<endl;
        continue;
      }
      int row, col;
      double maxDiff = (areas-areasGT).cwiseAbs().maxCoeff(&row, &col);
      if (maxDiff>tolerance){
        cout<<"The difference between areas and ground truth is "<<maxDiff<<" at ("<<row<<","<<col<<")"<<endl;
        continue;
      }
      maxDiff = (normals-normalsGT).cwiseAbs().maxCoeff(&row, &col);
      if (maxDiff>tolerance){
        cout<<"The difference between normals and ground truth is "<<maxDiff<<" at ("<<row<<","<<col<<")"<<endl;
        continue;
      }
      
      cout<<"Result is good!"<<endl;
      
    }
  }
}

