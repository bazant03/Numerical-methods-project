#include "GaussianElimination.h"

GaussianElimination::GaussianElimination(uint8 method, uint8 data):
    _method(method), _data(data) {
    if (data == 1) {
        this->assign_example_matrix(method);
        std::cout << "---------------------------------------------------------\n"
                  << "Przykładowa macierz:\n";
        this->print_matrix_A();
        std::cout << "---------------------------------------------------------\n";
    } else if (data == 2) {
        while (true) {
            std::cout << "---------------------------------------------------------\n"
                      << "Podaj wymiar macierzy kwadratowej (bez wektora wyrazów wolnych): \n"
                      << "---------------------------------------------------------\n";
            std::cin >> _dim;
            if (_dim > 0)
                break;
        }
        this->alloc_matrix_A();
        std::cout << "---------------------------------------------------------\n"
                  << "Podaj wartości elementów macierzy:\n";
        this->set_A();
        std::cout << "---------------------------------------------------------\n"
                  << "Podana macierz:\n";
        this->print_matrix_A();
        std::cout << "---------------------------------------------------------\n";
    }
}

GaussianElimination::~GaussianElimination() {
    if (_A != nullptr) {
        for (size_t i = 0; i < _dim; i++) {
            delete[] _A[i];
        }
        delete[] _A;
    }

    if (_x != nullptr) {
        delete[] _x;
    }

    if (_ind != nullptr) {
        delete[] _ind;
    }
}

void GaussianElimination::alloc_matrix_A() {
    _A = new double128*[_dim];
    for (size_t i = 0; i < _dim; i++) {
        _A[i] = new double128[_dim + 1];
    }
}

void GaussianElimination::assign_example_matrix(uint8 method) {
    switch (method) {
        case 1:
            _dim = 4;
            _A = new double128*[4]{
                            new double128[5]{1, 1, 0, -3, 1}, 
                            new double128[5]{1, 4, -1, -4, -2}, 
                            new double128[5]{0.5, 0.5, -3, -5.5, 1.5}, 
                            new double128[5]{1.5, 3, -5, -9, -0.5}
                        };
            break;
        case 2:
            _dim = 4;
            _A = new double128*[4]{
                            new double128[5]{2, 1, 1, 0, 4}, 
                            new double128[5]{4, 3, 3, 1, 6}, 
                            new double128[5]{8, 7, 9, 5, 8}, 
                            new double128[5]{6, 7, 9, 8, -2}
                        };
            break;
        case 3:
            /*_dim = 4;
            _A = new double128*[4]{
                            new double128[5]{2.25, -2.5, 4, -5.25, -1}, 
                            new double128[5]{-3, -7.5, 6.5, 0, 17}, 
                            new double128[5]{-6.25, -12.5, 0.25, 5.25, 24.25}, 
                            new double128[5]{9, 10, 7, -21, -33}
                        };*/
            _dim = 5;
            _A = new double128*[5]{
                            new double128[6]{14, -13, 3, -16, -42, -37}, 
                            new double128[6]{3.5, -18, 13, -23.75, -21, -5.5}, 
                            new double128[6]{3.5, 3, -5.25, 9.25, 10.5, 12.50}, 
                            new double128[6]{2, 14.5, -10.5, 18.5, 21, 23.5},
                            new double128[6]{1.5, 6.75, -9.25, 17, -10.5, -45.25}
                        };
            break;
    }
}

void GaussianElimination::set_A() {
    size_t i, j;
    for (i = 0; i < _dim; i++) {
        for (j = 0; j < _dim; j++) {
            std::cout << "a" << i + 1 << j + 1 << ": ";
            std::cin >> _A[i][j];
        }
        std::cout << "b" << i + 1 << ": ";
        std::cin >> _A[i][j];
    }
}

void GaussianElimination::print_matrix_A() {
    for (size_t i = 0; i < _dim; i++) {
        for (size_t j = 0; j < _dim + 1; j++) {
            printf("%7.2Lf", _A[i][j]);
        }
            std::cout << "\n";
    }
}

void GaussianElimination::alloc_vector_x() {
    _x = new double128[_dim];
}

void GaussianElimination::print_vector_x() {
    std::cout << "x = [";
    for (size_t i = 0; i < _dim; i++) {
        printf("%7.2Lf", _x[i]);
    }
    std::cout << "  ]T\n";
}

void GaussianElimination::set_vector_ind() {
    _ind = new size_t[_dim];
    for (size_t i = 0; i < _dim; i++) {
        _ind[i] = i;
    }
}

void GaussianElimination::swap_x() {
    double128 *x = new double128[_dim];
    for (size_t i = 0; i < _dim; i++) {
        x[_ind[i]] = _x[i];
    }
    delete [] _x;
    _x = x;
}

void GaussianElimination::reverse_substitution() {
    double128 sum;
    _x[_dim - 1] = _A[_dim - 1][_dim] / _A[_dim - 1][_dim - 1];
    for (size_t i = _dim - 1; i--;) {
        sum = 0.0;
        for (size_t j = i + 1; j < _dim; j++) {
            sum += _A[i][j] * _x[j];
        }
        _x[i] = (_A[i][_dim] - sum) / _A[i][i];
    }
}

void GaussianElimination::naive_elimination() {
    double128 r;
    for (size_t k = 0; k < _dim - 1; k++)  {
        
        for (size_t i = k + 1; i <= _dim - 1; i++) {
            
            if (fabs(_A[k][k]) <= _e) {
                _status = 1;
                _msg = "Na głównej przekątnej macierzy znajduje się 0,\nNie można rozwiązać układu metodą podstawową.\n";
                return;
            }

            r = _A[i][k] / _A[k][k];
            _A[i][k] = 0;
           
            for (size_t j = k + 1; j <= _dim; j++) {
                _A[i][j] -= r * _A[k][j];
            }

        }
    }
    _status = 0;
}

void GaussianElimination::partial_pivoting_elimination() {
    double128 r = 0, max = 0, temp = 0;
    size_t p;

    for (size_t k = 0; k < _dim - 1; k++)  {
        for (size_t i = k + 1; i <= _dim - 1; i++) {

            max = fabs(_A[k][k]);
            p = k;
            for (size_t l = k + 1; l < _dim; l++) {
                temp = fabs(_A[l][k]);
                if (temp > max) {
                    max = temp;
                    p = l; 
                }
            }

            if (max <= _e) {
                _status = 2;
                _msg = "Układ osobliwy";
                return;
            }

            for (size_t m = k; m < _dim + 1; m++) {
                    std::swap(_A[k][m], _A[p][m]);
            }

            r = _A[i][k] / _A[k][k];
            _A[i][k] = 0;
            for (size_t j = k + 1; j <= _dim; j++) {
                _A[i][j] -= r * _A[k][j];
            }
        }
    }
    _status = 0;
}

void GaussianElimination::complete_pivoting_elimination() {
    double128 r, max, temp;
    size_t p = 0, q = 0;

    for (size_t k = 0; k < _dim - 1; k++) {
        max = _e;

        for (size_t l = k; l < _dim; l++) {
            for (size_t m = k; m < _dim; m++) {
                temp = fabs(_A[l][m]);
                if (temp > max) {
                    max = temp;
                    p = l;
                    q = m;
                }
            }
        }

        if (max <= _e) {
            _status = 2;
            _msg = "Układ osobliwy";
            return;
        }

        for (size_t l = k; l < _dim + 1; l++) { 
            std::swap(_A[k][l], _A[p][l]);
        }
        
        for (size_t m = 0; m < _dim; m++) {
            std::swap(_A[m][k], _A[m][q]);
        }

        std::swap(_ind[k], _ind[q]);

        for (size_t i = k + 1; i < _dim; i++) {
            r = _A[i][k] / _A[k][k];
            _A[i][k] = 0;
            for (size_t j = k + 1; j < _dim + 1; j++) { 
                _A[i][j] -= r * _A[k][j];
            }
        }
    }
    _status = 0;
}
