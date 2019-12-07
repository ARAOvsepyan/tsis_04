#pragma once
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <numeric>
#include <vector>
#include <cmath>

using matrix = std::vector<std::vector<double>>;

const std::vector<std::string> trees = { "Oxford", "MSU", "MIPT", "TSU" };

const std::vector<std::string> criteria = { " Scholarship ", " Qualification ", " Cost of living ", " Prestigiousness " };

const std::vector<double> criteriaWeight = { 4, 6, 2, 8 };

std::vector<double> vectorNorm(const std::vector<double>& cw) {
    double sum = std::accumulate(cw.begin(), cw.end(), 0);

    std::vector<double> buff(cw.size());

    for (auto i = 0; i < (cw.size()); ++i)
        buff[i] = cw[i] / sum;

    return buff;
}

void printMatrix(const matrix& m) {
    std::cout << std::setw(10) << "";

    for (auto i = 0; i < criteria.size(); ++i) {
        std::cout << std::setw(20) << criteria[i];
    }

    std::cout << std::endl;

    for (auto i = 0; i < m.size(); ++i) {
        std::cout << std::setw(7) << trees[i];
        for (auto j = 0; j < m[i].size(); ++j) {
            std::cout << std::setw(19) << m[i][j];
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void printCritMatrix(const matrix& m) {
    std::cout << std::setw(15) << "";

    for (auto i = 0; i < m.size(); ++i) {
        std::cout << std::setw(24) << criteria[i];
    }

    std::cout << std::endl;

    for (auto i = 0; i < m.size(); ++i) {
        std::cout << std::setw(7) << criteria[i];

        for (auto j = 0; j < m[i].size(); ++j) {
            std::cout << std::setw(21) << m[i][j];
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void printCritPrioritet(const std::vector<double>& vec) {
    std::cout << "Criteria priority vector: ";

    for (auto i = 0; i < vec.size(); ++i) {
        std::cout << std::setw(7) << vec[i];
    }

    std::cout << std::endl << std::endl;
}

matrix matrixNorm(const matrix& m) {
    double sum = 0;
    std::vector<std::vector<double>> buff = m;

    for (auto i = 0; i < criteria.size(); ++i) {
        for (auto j = 0; j < m[i].size(); ++j) {
            sum += m[j][i];
        }
        for (auto j = 0; j < m[i].size(); ++j) {
            buff[j][i] = m[j][i] / sum; sum = 0;
        }
    }

    return buff;
}

void matrixNorm(matrix& m, size_t j) {
    double sum = 0;

    for (auto i = 0; i < 4; ++i) {
        sum += m[i][j];
    }

    for (auto i = 0; i < 4; ++i) {
        m[i][j] = m[i][j] / sum;
    }
}

void replaceCritConstr(const matrix& m, size_t k) {
    // Возьмем за главный критерий пристижность диплома

    std::vector<double> maxElements;

    for (auto j = 0; j < 4; ++j) {
        double max = m[0][j];
        for (auto i = 0; i < 4; ++i)
            if (m[i][j] > max)
                max = m[i][j];
        maxElements.push_back(max);
    }

    std::vector<double> minElements;

    for (auto j = 0; j < 4; ++j) {
        double min = m[0][j];
        for (auto i = 0; i < 4; ++i)
            if (m[i][j] < min)
                min = m[i][j];
        minElements.push_back(min);
    }

    auto normMatrix = m;

    for (auto i = 0; i < 4; ++i)
        for (auto j = 0; j < 4; ++j) {
            if (j == k) continue;
            normMatrix[i][j] = (normMatrix[i][j] - minElements[j]) / (maxElements[j] - minElements[j]);
        }

    maxElements.clear();

    for (auto j = 0; j < 4; ++j) {
        double max = normMatrix[0][j];
        for (auto i = 0; i < 4; ++i)
            if (normMatrix[i][j] > max)
                max = normMatrix[i][j];
        maxElements.push_back(max);
    }

    printMatrix(normMatrix);

    std::vector<double> solution;

    for (auto i = 0; i < 4; ++i) {
        double sum = 0; double max = normMatrix[0][i];
        if (normMatrix[i][0] >= 0.3 * maxElements[0] && normMatrix[i][1] >= 0.1 * maxElements[1] && normMatrix[i][2] >= 0.5 * maxElements[2]) {

            sum = std::accumulate(normMatrix[i].begin(), normMatrix[i].end(), 0);
            solution.push_back(sum);
        }
        else solution.push_back(0);
    }

    size_t num_items = std::count(solution.begin(), solution.end(), 0);

    if (num_items != solution.size()) {
        std::cout << "The value of the combined criterion: ";

        for (auto i = 0; i < solution.size(); ++i)
            std::cout << std::setw(7) << solution[i];

        std::cout << std::endl;
        auto max = std::max_element(solution.begin(), solution.end());
        std::cout << "The best solution: " << trees[std::distance(solution.begin(), max)] << std::endl;
    }
    else std::cout << "No solution found" << std::endl;
}

size_t manhattan(const std::pair<double, double>& a, const std::pair<double, double>& b) {
    return abs(a.first - b.first) + abs(a.second - b.second);
}

void metodPareto(const matrix& matrix, size_t k1, size_t k2) {
    std::vector<std::pair<double, double>> points;

    for (auto i = 0; i < 4; ++i)
        points.push_back(std::make_pair(matrix[i][k1], matrix[i][k2]));
    std::cout << "Points: ";
    for (const auto& it : points) std::cout << '(' << it.first << ", " << it.second << ") ";
    std::cout << std::endl;
    auto max_x = std::max_element(points.begin(), points.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.first < b.first;
    });
    auto max_y = std::max_element(points.begin(), points.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.second < b.second;
    });
    std::pair<double, double> utopiaPoint = std::make_pair(10, 10);
    std::cout << "Utopia points: (" << utopiaPoint.first << ", " << utopiaPoint.second << ") ";
    std::cout << std::endl; std::vector<double> solution;
    for (auto i = 0; i < points.size(); ++i) solution.push_back(manhattan(points[i], utopiaPoint));
    std::cout << "The value of the combined criterion: ";
    for (auto i = 0; i < solution.size(); ++i)
        std::cout << std::setw(7) << solution[i];
    std::cout << std::endl;
    auto min = std::min_element(solution.begin(), solution.end());

    std::cout << "The best solution: " << trees[std::distance(solution.begin(), min)] << std::endl;
}

void metodWeiComb(const matrix& m, const std::vector<double>& v) {
    auto buff = matrixNorm(m);
    printMatrix(buff);
    std::vector<double> solution;
    matrix critMatrix(m.size());

    for (auto i = 0; i < m.size(); ++i) {
        critMatrix[i].resize(m[i].size());
    }
    for (auto i = 0; i < m.size(); ++i) {
        for (auto j = 0; j < m[i].size(); ++j) {
            critMatrix[i][j] = v[i] / v[j];
        }
        std::cout << std::setw(70) << " Matrix for criterion" << std::endl;
        printCritMatrix(critMatrix);
    }

    std::vector<double> critPrioritet(critMatrix.size());
    for (auto i = 0; i < critMatrix.size(); ++i)
        for (auto j = 0; j < critMatrix[i].size(); ++j)
            critPrioritet[i] += critMatrix[i][j];
    critPrioritet = vectorNorm(critPrioritet);
    printCritPrioritet(critPrioritet);
    for (auto i = 0; i < m.size(); ++i) {
        double sum = 0; for (auto j = 0; j < m[i].size(); ++j) sum += buff[i][j] * critPrioritet[j];

        solution.push_back(sum);
    }
    std::cout << "The value of the combined criterion: ";
    for (auto i = 0; i < solution.size(); ++i)
        std::cout << std::setw(7) << solution[i];
    std::cout << std::endl;
    auto max = std::max_element(solution.begin(), solution.end());
    std::cout << "The best solution: " << trees[std::distance(solution.begin(), max)] << std::endl;
}

void printAlterVectorsPrioritet(const matrix& m) {
    std::cout << "Priority vectors of alternatives: " << std::endl;
    std::cout << std::endl;
    for (auto i = 0; i < m.size(); ++i) {
        std::cout << std::setw(7) << trees[i];
        for (auto j = 0; j < m[i].size(); ++j)
            std::cout << std::setw(7) << m[i][j];
        std::cout << std::endl;
    } std::cout << std::endl;
}

void printAlterMatrix(const matrix& m) {
    std::cout << std::setw(6) << "";
    for (auto i = 0; i < m.size(); ++i)
        std::cout << std::setw(10) << trees[i]; std::cout << std::endl;
    for (auto i = 0; i < m.size(); ++i) {
        std::cout << std::setw(7) << trees[i];
        for (auto j = 0; j < m[i].size(); ++j)
            std::cout << std::setw(10) << m[i][j];
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

matrix grtCritMatrix(const matrix& m, size_t k) {
    std::vector<double> cr;
    for (auto i = 0; i < m.size(); ++i)
        cr.push_back(m[i][k]);
    matrix critMatrix(m.size());
    for (auto i = 0; i < m.size(); ++i)
        critMatrix[i].resize(m[i].size());
    for (auto i = 0; i < m.size(); ++i)
        for (auto j = 0; j < m[i].size(); ++j)
            critMatrix[i][j] = cr[i] / cr[j];
    return critMatrix;
}

bool isAgreed(const matrix& m) {
    // оценка компонент собственного вектора
    std::vector<double> critVector;
    for (auto i = 0; i < m.size(); ++i) {
        double p = 1;
        for (auto j = 0; j < m[i].size(); ++j)
            p *= m[i][j];
        critVector.push_back(sqrt(sqrt(p)));
    }
    //Нормализация оценок
    double sum = std::accumulate(critVector.begin(), critVector.end(), 0);
    std::vector<double> norm_vec_cr(critVector.size());
    for (auto i = 0; i < (critVector.size()); ++i) norm_vec_cr[i] = critVector[i] / sum;
    //Вычисление максимального собственного числа матрицы
    std::vector<double> buffVector;
    for (int i = 0; i < 4; i++) {
        double sum = 0;
        for (int j = 0; j < 4; j++)
            sum += m[j][i];
        buffVector.push_back(sum * norm_vec_cr[i]);
    }
    auto max = std::accumulate(buffVector.begin(), buffVector.end(), 0);
    //Вычислим индекс согласования
    double is = (max - m.size()) / (m.size() - 1);
    //Оценка согласования для 4 значения случайной согласованности 0.9 
    double os = is / 0.9;
    std::cout << "assessment of consistency: " << os << std::endl;
    if (os <= 0.2)
        return true;
    return false;
}

void metodHierarchy(const matrix& m, const std::vector<double>& v) {
    std::vector<matrix> matrixs(m[0].size());
    std::vector<matrix> normMatrixs(m[0].size());
    for (auto i = 0; i < m[0].size(); ++i) {
        std::cout << " Matrix alternatives for criterion: " << criteria[i] << std::endl;
        matrixs[i] = grtCritMatrix(m, i);
        printAlterMatrix(matrixs[i]);
        isAgreed(matrixs[i]);
    }
    matrix prioritetAlters(m.size());
    for (auto k = 0; k < matrixs.size(); ++k) {
        std::vector<double> buff; for (auto i = 0; i < matrixs[k].size(); ++i) {
            double sum = 0;
            for (auto j = 0; j < matrixs[k][i].size(); ++j) {
                sum += matrixs[k][i][j];
            }
            buff.push_back(sum);
        }
        prioritetAlters[k] = vectorNorm(buff);
    }
    printAlterVectorsPrioritet(prioritetAlters);
    matrix critMatrix(m.size());
    for (auto i = 0; i < m.size(); ++i)
        critMatrix[i].resize(m[i].size());
    for (auto i = 0; i < m.size(); ++i)
        for (auto j = 0; j < m[i].size(); ++j)
            critMatrix[i][j] = v[i] / v[j];
    std::cout << std::setw(70) << " Matrix for criterion" << std::endl;
    printCritMatrix(critMatrix);
    std::vector<double> priority_criteria(critMatrix.size());
    for (auto i = 0; i < critMatrix.size(); ++i)
        for (auto j = 0; j < critMatrix[i].size(); ++j)
            priority_criteria[i] += critMatrix[i][j];
    priority_criteria = vectorNorm(priority_criteria);
    printCritPrioritet(priority_criteria);
    std::vector<double> solution;
    for (auto i = 0; i < prioritetAlters.size(); ++i) {
        double sum = 0;
        for (auto j = 0; j < priority_criteria.size(); ++j)
            sum += prioritetAlters[i][j] * priority_criteria[j];
        solution.push_back(sum);
    }
    std::cout << "The value of the combined criterion: ";
    for (auto i = 0; i < solution.size(); ++i)
        std::cout << std::setw(7) << solution[i];
    std::cout << std::endl;
    auto max = std::max_element(solution.begin(), solution.end());
    std::cout << "The best decision: " << trees[std::distance(solution.begin(), max)] << std::endl;
}