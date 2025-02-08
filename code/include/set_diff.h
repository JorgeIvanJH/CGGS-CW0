#ifndef SET_DIFF_HEADER_FILE
#define SET_DIFF_HEADER_FILE

#include <Eigen/Dense>
#include <set>

// Returns containing the elements of A that are not in B
Eigen::VectorXi set_diff(const Eigen::VectorXi& A,
                         const Eigen::VectorXi& B) {

    //By converting B into a set, we can quickly check whether an element of A is in B
    std::set<int> setB;
    for (int i = 0; i < B.size(); ++i) {
        setB.insert(B[i]);
    }
    
    //If the element is in A, but not in B, add it to the result vector
    std::vector<int> result;
    for (int i = 0; i < A.size(); ++i) {
        if (setB.find(A[i]) == setB.end()) {
            result.push_back(A[i]);
        }
    }
    
    // Convert result to Eigen::VectorXi
    Eigen::VectorXi diff(result.size());
    for (size_t i = 0; i < result.size(); ++i) {
        diff[i] = result[i];
    }
    
    return diff;
}

#endif
