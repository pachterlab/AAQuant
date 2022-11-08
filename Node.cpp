#include "Node.hpp"

Node::Node() {
    read_counts = 0;
    abundance.clear();
}

void Node::clear(const UnitigMap<Node> &um_dest) {
    read_counts = 0;
    abundance.clear();
    r = Roaring();
}

// TODO:
// Take a better look at this -- might be setting counts to zero in some cases
void Node::concat(const UnitigMap<Node> &um_dest, const UnitigMap<Node> &um_src) {
    Node *data_dest = um_dest.getData();
    Node *data_src = um_src.getData();

    r = Roaring();

    float w_dest = um_dest.len / (um_dest.len + um_src.len);
    float w_src = um_src.len / (um_dest.len + um_src.len);

    // Fill abundance vector with zeros
    std::vector<std::pair<uint32_t, float> > elems;
    data_dest->abundance.getElements(elems);
    abundance = SparseVector<float>();
    for (const auto &it : elems) {
        abundance.insert(it.first, it.second*w_dest);
    }
    elems.clear();
    data_src->abundance.getElements(elems);
    for (const auto &it : elems) {
        if (abundance.contains(it.first)) {
            abundance[it.first] += it.second*w_src;
        } else {
            abundance.insert(it.first, it.second*w_src);
        }
    }
    /*
    abundance.resize(data_dest->abundance.size(), 0.);

    for (size_t idx = 0; idx < data_dest->abundance.size(); ++idx) {
        abundance[idx] += data_dest->abundance[idx]*w_dest + data_src->abundance[idx]*w_src;
    }
    */
    read_counts = data_dest->read_counts*w_dest + data_src->read_counts*w_src;

    // Concatenate deletion flag bit vectors.
    // r = [um_dest.r, um_src.r] if (um_dest.strand && um_src.strand)
    // r = [rev(um_dest.r), um_src.r] if (!um_dest.strand && um_src.strand)
    // r = [um_dest.r, rev(um_src.r)] if (um_dest.strand && !um_src.strand)
    // r = [rev(um_dest.r), rev(um_src.r)] if (!um_dest.strand && !um_src.strand) 

    if (!data_dest->r.isEmpty()) {
        if (um_dest.strand) {
            // Forward strand copy
            r.addRange(0, um_dest.len);
        } else {
            // Reverse strand copy
            for (int i = um_dest.len-1; i >= 0; --i) {
                // Check whether bit is set
                if (data_dest->r.contains(um_dest.len - 1 - i)) {
                    r.add(i);
                }
            }
        }
    }
 
    if (!data_src->r.isEmpty()) {
        if (um_src.strand) {
            // Forward strand shift and copy
            for (int i = um_dest.len; i < um_dest.len + um_src.len; ++i) {
                if (data_src->r.contains(i - um_dest.len)) {
                    r.add(i);
                }
            }
        } else {
            // Reverse strand shift and copy
            for (int i = um_dest.len + um_src.len - 1; i >= um_dest.len; --i) {
                if (data_src->r.contains(um_dest.len + um_src.len - 1 - i)) {
                    r.add(i);
                }
            }
        }
    }
} 

void Node::extract(const UnitigMap<Node> &um_src, bool last_extraction) {
    Node *data = um_src.getData();

    abundance = SparseVector<float>(data->abundance);
    read_counts = data->read_counts;
    // Copy delete bit vetor
    if (!data->r.isEmpty()) {
        for (size_t i = 0; i < um_src.len; ++i) {
            if (data->r.contains(um_src.dist + i)) {
                r.add(i);
            }
        }
    }
}

std::string Node::toString() const {
    return "Total abundance: " + std::to_string(read_counts);
}

// Flags the kmers from x (inclusive) to y (noninclusive)
void Node::flag(const uint64_t x, const uint64_t y) {
    r.addRange(x, y);
}

void Node::flag(const UnitigMap<Node> &um) {
    r.addRange(0, um.len);
}
