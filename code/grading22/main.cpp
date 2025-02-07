#include "readOFF.h"
#include "create_edge_list.h"
#include "serialization.h"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <filesystem>

using namespace Eigen;
using namespace std;

//The mesh quantities
MatrixXi F, E, H;
MatrixXd V;
VectorXi boundEMask;

namespace fs = std::filesystem;



int main()
{
  
  std::string folderPath(DATA_PATH); // Replace with your folder path
  for (const auto& entry : fs::directory_iterator(folderPath)) {
    if (entry.is_regular_file() && entry.path().extension() == ".off") {
      cout<<"Working on file "<<entry.path().filename()<<endl;
      std::string dataName = entry.path().string();
      dataName.erase(dataName.size() - 4, 4);
      std::ifstream ifs(dataName+"-section22.data", std::ofstream::binary);
      MatrixXi HGT, EGT;
      VectorXi boundEMaskGT;
      readOFF(entry.path().string(), V, F);
      create_edge_list(F,H, E, boundEMask);
      deserializeMatrix(HGT, ifs);
      deserializeMatrix(EGT, ifs);
      deserializeVector(boundEMaskGT, ifs);
      ifs.close();
      
      //comparisons and feedback
      if ((H.rows()!=HGT.rows())||(H.cols()!=HGT.cols())){
        cout<<"H is of the wrong size!"<<endl;
        continue;
      }
      if ((E.rows()!=EGT.rows())||(E.cols()!=EGT.cols())){
        cout<<"E is of the wrong size!"<<endl;
        continue;
      }
      if ((boundEMask.rows()!=boundEMaskGT.rows())||(boundEMask.cols()!=boundEMaskGT.cols())){
        cout<<"boundEMask is of the wrong size!"<<endl;
        continue;
      }
      int row, col;
      int maxDiff = (H-HGT).cwiseAbs().maxCoeff(&row, &col);
      if (maxDiff>0){
        cout<<"The difference between H and ground truth is "<<maxDiff<<" at ("<<row<<","<<col<<")"<<endl;
        continue;
      }
      maxDiff = (E-EGT).cwiseAbs().maxCoeff(&row, &col);
      if (maxDiff>0){
        cout<<"The difference between E and ground truth is "<<maxDiff<<" at ("<<row<<","<<col<<")"<<endl;
        continue;
      }
      
      maxDiff = (boundEMask-boundEMaskGT).cwiseAbs().maxCoeff(&row, &col);
      if (maxDiff>0){
        cout<<"The difference between boundEMask and ground truth is "<<maxDiff<<" at ("<<row<<","<<col<<")"<<endl;
        continue;
      }
      
      
      cout<<"Result is good!"<<endl;
      
    }
  }
  
}

