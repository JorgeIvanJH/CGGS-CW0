#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include "readOFF.h"
#include "compute_areas_normals.h"
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;
using namespace polyscope;

//The mesh quantities
MatrixXi F;
MatrixXd V;


int main()
{
    cout << "OK" << endl;
    readOFF(DATA_PATH "/Kitten.off",V, F);
    polyscope::init();
    SurfaceMesh* pSurfaceMesh = polyscope::registerSurfaceMesh("mesh", V, F);
    VectorXd areas;
    MatrixXd normals;
    compute_areas_normals(V,F, areas, normals);
    pSurfaceMesh->addFaceVectorQuantity("face normals", normals);
    pSurfaceMesh->addFaceScalarQuantity("face areas", areas);
    polyscope::show();
}

