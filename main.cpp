#include "GaussianElimination.h"

using namespace std;

uint8 choose_data() {
    uint16 t;
    while (true) {
        cout << "Wybierz rodzaj danych: \n"
            << "1.Przykładowe dane z programu\n"
            << "2.Wprowadź własne dane\n"
            << "---------------------------------------------------------\n";
        
        cin >> t;

        if (t == 1 || t == 2)
            break;
    }
    return static_cast<uint8>(t);
}

uint8 choose_method() {
    uint16 t;
    while (true) {
        cout << "Wybierz metodę do rozwiązania układu równań: \n"
         << "1.metoda podstawowa\n"
         << "2.metoda z częściowym wyborem elementu maksymalnego\n"
         << "3.metoda z pełnym wyborem elementu maksymalnego\n"
         << "---------------------------------------------------------\n";
    
        cin >> t;

        if (t >= 1 && t <= 3)
            break;
    }
    return static_cast<uint8>(t);
}

int main() {

    std::cout << "---------------------------------------------------------\n";
    uint8 method = choose_method();
    std::cout << "---------------------------------------------------------\n";
    uint8 data = choose_data();

    GaussianElimination *g = new GaussianElimination(method, data);

    switch (g->_method) {
        case 1:
            g->naive_elimination();
            break;
        case 2:
            g->partial_pivoting_elimination();
            break;
        case 3:
            g->set_vector_ind();
            g->complete_pivoting_elimination();
            break;
    }

    if (g->_status == -1) {
        cout << "Nie udało się rozwiązać układu\n";
             delete g;
        return 0;
    } else if (g->_status > 0) {
        cout << g->_msg << "\n";
             delete g;
        return 0;
    }

    g->alloc_vector_x();
    g->reverse_substitution();

    if (g->_method == 3) { 
        g->swap_x();
    }
    
    std::cout << "Rozwiązanie:\n";   // P1 x = [-1 -1 1 -1]T
    g->print_vector_x();             // P2 x = [2 -1 1 -2] T
                                     // P3.1 x = [-1 -1 1 1] T
    delete g;

    return 0;
}