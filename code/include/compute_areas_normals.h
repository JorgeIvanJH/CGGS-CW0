#ifndef COMPUTE_AREAS_NORMALS_HEADER_FILE
#define COMPUTE_AREAS_NORMALS_HEADER_FILE

#include <Eigen/Dense>

void compute_areas_normals(const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    Eigen::VectorXd& areas,
    Eigen::MatrixXd& normals) {
    using namespace Eigen;
    using namespace std;
    areas.resize(F.rows()); // areas set to #rows in F and 1 column
    normals.resize(F.rows(), 3); // normals set to #rows in F and 3 columns
    cout << "V. shape: " << V.rows() << " " << V.cols() << endl;
    cout << "F. shape: " << F.rows() << " " << F.cols() << endl;


    //TODO

    // LOOP APPROACH
    for (int i = 0; i < F.rows(); i++) {
        RowVector3d v0 = V.row(F(i, 0));
        RowVector3d v1 = V.row(F(i, 1));
        RowVector3d v2 = V.row(F(i, 2));
        RowVector3d e0 = v1 - v0;
        RowVector3d e1 = v2 - v0;
        Vector3d nf = e0.cross(e1);
        double area = 1 / 2 * nf.norm();
        areas(i) = area;
        normals.row(i) = nf.normalized();
    }


    // MAP APPROACH
    // MatrixXd e0(F.rows(), 3);
    // MatrixXd e1(F.rows(), 3);
    // cout << "OK 0" << endl;
    // for (int i = 0; i < F.rows(); i++) {
    //     e0.row(i) = V.row(F(i, 1)) - V.row(F(i, 0));
    //     e1.row(i) = V.row(F(i, 2)) - V.row(F(i, 0));
    // }
    // normals = e0.rowwise().cross(e1);
    // areas = 0.5 * normals.rowwise().norm();
    // normals = normals.rowwise().normalized();



    // cout << "areas: " << areas << endl;
    // cout << "normals: " << normals << endl;
}


#endif
