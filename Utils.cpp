#include "Utils.hpp"

void H1_neighborhood(const std::string &unitig, std::vector<std::string> &h1) {
    char bases[4] = {'A', 'C', 'G', 'T'};
    h1.reserve(unitig.size()*3);
    std::string ham;
    for (size_t i = 0; i < unitig.length(); ++i) {
        ham = unitig;
        for (auto& c : bases) {
            if (c != unitig[i]) {
                ham[i] = c;
                h1.push_back(ham);
            }
        }
    }
}
