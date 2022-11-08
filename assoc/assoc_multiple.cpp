#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <cmath>

#define REMORA_USE_SIMD

#include "CLI11.hpp"
#include "remora/include/remora/remora.hpp"


void read_key_value(const std::string path, std::vector<std::string> &k, std::vector<float> &v) {
    std::cout << "Reading markers from " << path << std::endl;
    std::string pn;
    float val;
    std::ifstream infile(path);
    while (infile >> pn >> val) {
        k.push_back(pn);
        v.push_back(val);
    }
    infile.close();
}

void parse_coord_row_col(std::string &path, std::vector<std::string> &pns, std::vector<std::string> &kmers) {
    std::set<std::string> _pns, _kmers;
    std::string pn, kmer;
    float count;
    std::ifstream infile(path);
    while (infile >> pn >> kmer >> count) {
        _pns.insert(pn);
        _kmers.insert(kmer);
    }
    pns.assign(_pns.begin(), _pns.end());
    kmers.assign(_kmers.begin(), _kmers.end());
}

void square_inverse(std::vector<float> &vals, remora::matrix<float> &XTXI) {
    float sum = 0, sumsq = 0, n = vals.size();

    for (const auto &val : vals) {
        sum += val;
        sumsq += (val*val);
    }

    float denom = 1 / ((n * sumsq) - (sum * sum));
    XTXI(0, 0) = sumsq * denom;
    XTXI(0, 1) = -(sum * denom);
    XTXI(1, 0) = XTXI(0, 1);
    XTXI(1, 1) = n * denom;
    std::cout << "XTXI:" << std::endl;
    std::cout << XTXI(0, 0) << "\t" << XTXI(0, 1) << std::endl << XTXI(1, 0) << "\t" << XTXI(1, 1) << std::endl;
}

void read_coord_mtx(remora::matrix<float> &Y, std::string &path,
                    std::unordered_map<std::string, size_t> &pn2idx,
                    std::unordered_map<std::string, size_t> &kmer2idx) {
    std::string pn, kmer;
    float count;
    std::ifstream infile(path);
    while (infile >> pn >> kmer >> count) {
        if (pn2idx.contains(pn) && kmer2idx.contains(kmer)) {
            Y(pn2idx[pn], kmer2idx[kmer]) = count;
        }
    }
    infile.close();
}

void t_vals(remora::matrix<float> &betas,
            remora::vector<float> &se,
            remora::vector<float> &t) {
    for (int i = 0; i < se.size(); ++i) {
        t(i) = std::abs(betas(1, i)) / se(i);
    }
}

void apply_threshold(remora::vector<float> &t,
                     float t_threshold,
                     std::vector<int> &indices) {
    indices.clear();
    for (int i = 0; i < t.size(); ++i) {
        if (t[i] > t_threshold) {
            indices.push_back(i);
        }
    }
}

void std_err(remora::matrix<float> &X,
             remora::matrix<float> &Y,
             remora::matrix<float> &betas,
             remora::vector<float> &se) {

    size_t n = X.size1();
    float x_bar = remora::sum(column(X, 1)) / n;

    float denom = 0;
    for (size_t i = 0; i < n; ++i) {
        float sub = X(i, 1) - x_bar;
        denom += sub*sub;
    }
    denom = std::sqrt(denom);

    float num;
    for (size_t i = 0; i < Y.size2(); ++i) {
        num = 0;
        for (size_t j = 0; j < n; ++j) {
            // (y - y_hat)^2
            float sub = Y(j, i) - (betas(1, i)*X(j, 1) + betas(0, i));
            num += sub*sub;
        }
        num = std::sqrt(num / (n-2));
        se(i) = num / denom;
    }
}

void write(std::string &path,
           std::vector<std::string> &y_kmers,
           remora::matrix<float> &betas,
           remora::vector<float> &se,
           int df,
           std::vector<int> &indices) {

    std::cout << "Writing results to " << path << std::endl;

    std::ofstream outfile(path);
    outfile << "kmer\tcoef\tstderr\tdf" << std::endl;
    //for (size_t i = 0; i < y_kmers.size(); ++i) {
    for (const auto i : indices) {
        outfile << y_kmers[i] << "\t";
        outfile << betas(1, i) << "\t";
        outfile << se(i) << "\t";
        outfile << df << std::endl;
    }
    outfile.close();
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    // ============================
    std::string target_path, counts_path, out_path;
    float threshold = -1.;
    CLI::App app{"Run associations for unitig abundances, with an optional t-value threhsold."};
    app.add_option("-x", target_path, "Association target.")->required();
    app.add_option("-y,", counts_path, "Unitig abundances.")->required();
    app.add_option("-o,--out", out_path, "Path to output directory.")->required();
    // Threshold (R): qt((1 - (0.05/20000)/2), 13000)
    app.add_option("-t,--threshold", threshold, "T-value threshold. Only write out associations with a T-value higher than the threshold.");
    CLI11_PARSE(app, argc, argv);

    // Check whether target_path contains multiple markers or just a single one
    // ========================================================================
    std::vector<std::string> x_paths;
    std::ifstream infile(target_path);
    std::string line;
    std::getline(infile, line);
    if (line.substr(0, 1) == "/") {
        x_paths.push_back(line);
        while (!infile.eof()) {
            std::getline(infile, line);
            x_paths.push_back(line);
        }
    } else {
        x_paths.push_back(target_path);
    }
    infile.close();

    // Read Y column and row names
    // ===========================
    std::vector<std::string> y_pns, y_kmers;
    std::unordered_map<std::string, size_t> pn2idx, kmer2idx;

    parse_coord_row_col(counts_path, y_pns, y_kmers);

    // Set up Y
    // ========
    pn2idx.reserve(y_pns.size());
    kmer2idx.reserve(y_kmers.size());
    for (int i = 0; i < y_pns.size(); ++i) pn2idx[y_pns[i]] = i;
    for (int i = 0; i < y_kmers.size(); ++i) kmer2idx[y_kmers[i]] = i;

    remora::matrix<float> Y(y_pns.size(), y_kmers.size(), 0.);
    read_coord_mtx(Y, counts_path, pn2idx, kmer2idx);


    std::vector<std::string> pns;
    std::vector<float> t_val;

    for (const auto path : x_paths) {
        pns.clear();
        t_val.clear();
        // Read X data
        // ===========
        read_key_value(path, pns, t_val);

        // Set up X and (X^TX)^-1
        // ======================
        remora::matrix<float> XTXI(2, 2);
        // Have X be in R^{pns.size(), 1} for normalized phenotypes
        // Remember to rank transform the subset of phenotypes to be associated.
        remora::matrix<float> X(pns.size(), 2, 1.);

        square_inverse(t_val, XTXI);

        for (int i = 0; i < t_val.size(); ++i) {
            X(i, 1) = t_val[i];
        }

        remora::matrix<float> betas = (XTXI % trans(X)) % Y;

        remora::vector<float> se(Y.size2());
        std_err(X, Y, betas, se);

        // Threshold on t-values
        remora::vector<float> t(se.size());
        t_vals(betas, se, t);
        std::vector<int> indices;
        apply_threshold(t, threshold, indices);

        int df = X.size1() - 2;

        std::string file_name = out_path + path.substr(path.find_last_of("/"));
        write(file_name, y_kmers, betas, se, df, indices);
    }
}
