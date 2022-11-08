#include <iostream>

#include "SparseVector.hpp"

void print_vector(SparseVector<float> &sv, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << i << ":" << sv.get(i, 0);
        if (i < size-1) std::cout << ", ";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    SparseVector<float> sv;
    std::vector<std::pair<uint32_t, float> > elems;
    sv.getElements(elems);
    for (const auto &it : elems) {
        std::cout << it.first << ": " << it.second << std::endl;
    }
    elems.clear();
    for (int i = 0; i < 10; i++) {
        sv.insert(i*2, i);
    }
    print_vector(sv, 20);
    sv.getElements(elems);
    for (const auto &it : elems) {
        std::cout << it.first << ": " << it.second << std::endl;
    }
    elems.clear();

    sv.insert(11, 400.0);
    sv.insert(17, 5000.0);
    print_vector(sv, 20);
    sv.getElements(elems);
    for (const auto &it : elems) {
        std::cout << it.first << ": " << it.second << std::endl;
    }
    elems.clear();

    sv.remove(11);
    print_vector(sv, 20);

    sv.getElements(elems);
    for (const auto &it : elems) {
        std::cout << it.first << ": " << it.second << std::endl;
    }

    std::cout << std::endl;
}
