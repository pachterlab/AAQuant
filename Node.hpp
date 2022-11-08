#ifndef ACDBG_NODE_HPP
#define ACDBG_NODE_HPP
#include <vector>

#include <bifrost/CompactedDBG.hpp>

#include "SparseVector.hpp"

class Node: public CDBG_Data_t<Node> {
    public:
        Node();
        void clear(const UnitigMap<Node> &um_dest);
        void concat(const UnitigMap<Node> &um_dest, const UnitigMap<Node> &um_src);
        void extract(const UnitigMap<Node> &um_src, bool last_extraction);
        std::string toString() const;

        void flag(const uint64_t x, const uint64_t y);
        void flag(const UnitigMap<Node> &um);

        SparseVector<float> abundance;
        float read_counts;
        // bitmask to flag individual constituent kmers for deletion
        Roaring r;
};
#endif
