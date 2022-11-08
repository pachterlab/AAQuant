#include "RemoveTips.hpp"
/******************************************************
 *                                                    *
 *   High cardinality tips removal helper functions   *
 *                                                    *
 ******************************************************/
// Returns
//      number of predecessors, if um has predecessors and no successors,
//      number of successors, if um has successors and no predecessors,
//      0, otherwise.
int is_partial_tip_group(const UnitigMap<Node> &um) {
    size_t n_succ = um.getSuccessors().cardinality();
    size_t n_pred = um.getPredecessors().cardinality();
    if (n_pred != 0 && n_succ == 0) {
        return n_pred;
    } else if (n_pred == 0 && n_succ != 0) {
        return n_succ;
    }
    return 0;
}

// Returns
//      cardinality, if um has cardinality predecessors and no successors,
//      cardinality, if um has cardinality successors and no predecessors,
//      0, otherwise.
int is_partial_tip_group(const UnitigMap<Node> &um, int cardinality) {
    size_t n_succ = um.getSuccessors().cardinality();
    size_t n_pred = um.getPredecessors().cardinality();
    if ((!n_pred || !n_succ) && (n_pred + n_succ == cardinality)) {
        return cardinality;
    }
    return 0;
}

bool has_tip_neighbors(const UnitigMap<Node> &um, int cardinality) {
    for (auto& succ : um.getSuccessors()) {
        if (is_partial_tip_group(succ)) return 1;
    }
    for (auto& pred : um.getPredecessors()) {
        if (is_partial_tip_group(pred)) return 1;
    }
    return 0;
}

bool is_tip(const UnitigMap<Node> &um) {
    int c = is_partial_tip_group(um);
    if (c) {
        return (um.mappedSequenceToString().length() == 31) || has_tip_neighbors(um, c);
    }
    return 0;
}
