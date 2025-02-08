#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include "readOFF.h"
#include "unique.h"
#include "sort_rows.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <set>
#include <array>
#include "create_edge_list.h"
#include "slice_columns_sparse.h"
#include "set_diff.h"
#include "harmonic_interpolation.h"

using namespace Eigen;
using namespace std;

//The mesh quantities
MatrixXi F;
MatrixXd V;
MatrixXi E, H;
VectorXi boundEMask;

int main()
{
	cout << "BEGIN - Bunny -" << endl;
    readOFF(DATA_PATH "/Bunny.off",V, F);
    polyscope::init();
    polyscope::SurfaceMesh* pSurf = polyscope::registerSurfaceMesh("mesh", V, F);
    create_edge_list(F, H, E, boundEMask);
    
    VectorXi B(2); B<<1000, 7642; // subset of vertex indices that represents the fixed vertices in the interpolation problem
    VectorXd xB(2); xB<<0.0, 1.0; // values of the fixed vertices
    VectorXd x = harmonic_interpolation(V, E, B, xB);
    
    //Visualization
    MatrixXd VB(B.size(),3);
    for (int i=0;i<B.size();i++)
        VB.row(i) = V.row(B(i));
    
    pSurf->addVertexScalarQuantity("harmonic function", x);
    polyscope::registerPointCloud("Fixed Points", VB);

    polyscope::show();
    
}

