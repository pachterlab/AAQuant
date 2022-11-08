#ifndef REMOVE_TIPS_HPP
#define REMOVE_TIPS_HPP
#include <bifrost/UnitigMap.hpp>

#include "Node.hpp"

int is_partial_tip_group(const UnitigMap<Node>& um);
int is_partial_tip_group(const UnitigMap<Node>& um, int cardinality);
bool has_tip_neighbors(const UnitigMap<Node>& um, int cardinality);
bool is_tip(const UnitigMap<Node>& um);
#endif
