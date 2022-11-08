#ifndef SPARSE_VECTOR_HPP
#define SPARSE_VECTOR_HPP
#include <vector>
#include <bifrost/roaring.hh>

template <class T>
class SparseVector {
    public:
        SparseVector();
        SparseVector(SparseVector<T> &&arg);
        SparseVector(const SparseVector<T> &arg);
        SparseVector<T>& operator=(SparseVector<T> &&other);
        SparseVector<T>& operator=(const SparseVector<T> &other);
        void insert(size_t i, const T &elem);
        void insert(const std::pair<size_t, T> &elem);
        T remove(size_t i);
        void clear();
        // Populates elems with (index, element) tuples
        void getElements(std::vector<std::pair<uint32_t, T> > &elems) const;
        T get(size_t i) const;
        T get(size_t i, const T& def) const;
        bool contains(size_t i) const;
        bool isEmpty() const;

        T& operator[] (size_t i);
        const T& operator[] (size_t i) const;

        // BEGIN DEBUG
        std::vector<T> v;
        // END DEBUG

    private:
        Roaring r;
        //std::vector<T> v;
};

#include "SparseVector.tcc"
#endif
