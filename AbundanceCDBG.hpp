#ifndef ABUNDANCE_CDBG_HPP
#define ABUNDANCE_CDBG_HPP

#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <random>

#include <bifrost/CompactedDBG.hpp>

//#include "BAMReader.hpp"
#include "Utils.hpp"
#include "RemoveTips.hpp"


struct ACDBG_Build_opt {
    bool verbose;
    bool debug;
    size_t nb_threads;

    size_t count_threshold;
    float ratio_threshold;

    ACDBG_Build_opt() : verbose(false), debug(false), nb_threads(1), count_threshold(50),
                        ratio_threshold(.25) {}
};

class AbundanceCDBG {
    private:
        CompactedDBG<Node> dbg;

    public:
        std::vector<std::string> pns;
        std::unordered_map<std::string, size_t> pn2idx;

        AbundanceCDBG();

        bool read(const std::string &dbg_path, bool verbose);
        bool build(const ACDBG_Build_opt &opt, const std::vector<std::string> &bam_paths, const std::string &region, bool verbose);
        bool build(const ACDBG_Build_opt &opt, const std::string &fasta_path, bool verbose);
        bool write(const std::string path, size_t nb_threads, bool verbose);

        //void resize_abundance_vectors(size_t size);
        // remove_kmers takes a vector of kmers and removes them from the graph
        // while subtract_kmers takes a path to a file containing one kmer per
        // line
        bool remove_kmers(std::vector<std::string> &kmers_to_remove, bool verbose);
        bool remove_kmers(std::vector<Kmer> &kmers_to_remove, bool verbose);
        bool subtract_graph(const std::string &path, size_t threads, bool verbose);
        bool subtract_kmers(const std::string &path, bool verbose);
        void detect_tips(std::vector<Kmer> &tips);
        void remove_tips(bool verbose);

        void add_noise(float size);

        bool prune_individual_mean(const std::string &path, bool flag, std::vector<Kmer> &canary, bool keep, bool verbose);
        bool remove_low_wrt_h1(size_t count_threshold, float ratio_threshold, bool flag, std::vector<Kmer> &canary, bool keep, bool verbose);

        void flag_kmers(std::vector<Kmer> &kmers_to_flag, bool verbose);
        void flag_unitigs(std::vector<Kmer> &unitigs_to_flag, bool verbose);
        void flag_remove_tips(bool verbose);
        bool flag_subtract_graph(const std::string &path, size_t threads, bool verbose);
        bool flag_subtract_kmers(const std::string &path, bool verbose);
        bool reload_unflagged(ACDBG_Build_opt &opt, std::string od, bool verbose);
        bool remove_flagged_unitig(UnitigMap<Node> &um, std::vector<Kmer> &canary, bool keep, bool verbose);
        bool remove_flagged_unitigs(std::vector<Kmer> &canary, bool keep, bool verbose);

        // Not normalized
        void write_abundances(std::string path);
        // Normalized
        void write_abundances(std::string path, std::vector<float> &idx2total);
        void write_kmer2unitig(std::string path);

        void count_kmers(const std::string &fasta_path, size_t len_pn, bool verbose);
        void count_kmers(const std::unordered_map<std::string, std::string> &bam_files, const std::string &region, bool verbose);
        void avg_counts(std::vector<Kmer> &canary, bool keep, bool flag, bool verbose);

        bool simplify(bool verbose);

        // Total number of counts per pn
        void tally_counts(std::vector<float> &idx2total);

        // Sanity checks
        // Allows you to check the retention of a specified set of kmers after
        // the various pruning steps
        float canary_kmer_retention(std::vector<Kmer> &kmers);
        bool is_in_kmer_list(const UnitigMap<Node> &um, std::vector<Kmer> &kmers);
};
#endif
