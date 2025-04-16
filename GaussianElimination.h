#ifndef GAUSSIANELIMINATION_H
#define GAUSSIANELIMINATION_H

#include <cmath>
#include <iostream>

typedef long double double128;
typedef unsigned char uint8;
typedef unsigned short int uint16;

class GaussianElimination {
    public:
        size_t _dim = 0;          // Wymiar macierzy 1..n
        double128 **_A = nullptr; 
        double128 *_x = nullptr;   
        size_t *_ind = nullptr;   // Wektor indeksów 
        int _status = -1;
        uint8 _method = 0;
        uint8 _data = 0;
        const double128 _e = 0.000'000'1;
        std::string _msg;

        GaussianElimination(uint8 method, uint8 data);
        ~GaussianElimination();  

        void alloc_matrix_A();
        void assign_example_matrix(uint8 method);
        void set_A();
        void print_matrix_A();
        void alloc_vector_x();
        void print_vector_x();
        void set_vector_ind();
        void print_vectior_ind();
        void swap_x();

        void reverse_substitution();
        void naive_elimination();
        void partial_pivoting_elimination();
        void complete_pivoting_elimination();
};

#endif
