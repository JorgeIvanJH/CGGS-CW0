#ifndef CREATE_EDGE_LIST_HEADER_FILE
#define CREATE_EDGE_LIST_HEADER_FILE

#include <Eigen/Dense>
#include "unique.h"
#include "sort_rows.h"

using namespace std;
using namespace Eigen;

void create_edge_list(const Eigen::MatrixXi& F,
                      Eigen::MatrixXi& H,
                      Eigen::MatrixXi& E,
                      Eigen::VectorXi& boundEMask){
    int num_faces = 5;//F.rows();
  H.resize(3* num_faces, 2);  //halfedges
  E.resize(0,2);
  boundEMask.resize(0);
  //TODO
  // halfedges (H): size 3|F| × 2
  // edges (E): size |E| × 2
  for (int i = 0; i < num_faces; i++)
  {
    for (int j = 0; j < 3; j++)
    {
		H.row(3 * i + j) << F(i, j), F(i, (j + 1) % 3);
      
    }
  }
  Eigen::MatrixXi H_copy = H;
  sort_rows(H_copy);



  cout << "H: " << H << endl;
  cout << "H_copy: " << H_copy << endl;

   vector<int> uniqueIndices, counts, inverseIndices;

   unique(H_copy,uniqueIndices, counts, inverseIndices);

   for (int i = 0; i < uniqueIndices.size(); i++)
   {
	   cout << "uniqueIndices: " << uniqueIndices[i] << endl;
   }

  // E.resize(uniqueIndices.size(), 2);
  // boundEMask.resize(H.rows(), 1);
  // for (int i = 0; i < uniqueIndices.size(); i++)
  // {
  //   E.row(i) = H_copy.row(uniqueIndices[i]);
  //   boundEMask(uniqueIndices[i]) = counts[i] == 1;
  // }


}



#endif
