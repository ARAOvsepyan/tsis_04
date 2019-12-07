#include "lab_04.h"

int main() {
    std::cout << std::fixed << std::setprecision(3) << std::right;
    matrix m = { { 1, 10, 2, 10 },
                 { 7, 8, 4, 9 },
                 { 4, 5, 6, 7 },
                 { 3, 1, 10, 1 }
    };
    printMatrix(m);
    std::cout << "#######################################################################################" << std::endl; std::cout << std::endl;
    std::cout << std::setw(60) << "1. Replacing criteria with constraints" << std::endl << std::endl;
    replaceCritConstr(m, 3);
    std::cout << "#######################################################################################" << std::endl;
    std::cout << std::endl;
    std::cout << std::setw(63) << "2. Formation and judgment of the Pareto set" << std::endl << std::endl;
    metodPareto(m, 1, 2);
    std::cout << "#######################################################################################" << std::endl;
    std::cout << std::endl;
    std::cout << std::setw(62) << "3. Criteria weighing and combining method" << std::endl << std::endl;
    metodWeiComb(m, criteriaWeight);
    std::cout << "#######################################################################################" << std::endl;
    std::cout << std::endl;
    std::cout << std::setw(57) << "4. Hierarchy analysis method" << std::endl << std::endl;
    metodHierarchy(m, criteriaWeight);
    std::cout << "#######################################################################################" << std::endl;
    return 0;
}