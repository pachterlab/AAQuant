#include <string>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <cstdlib>

#include "CLI11.hpp"

#include "AbundanceCDBG.hpp"

const int LEN_PN=7;

int main(int argc, char* argv[]) {
    // Argument declarations
    bool verbose = false;
    bool dry_run = false;
    int threads = 1;
    std::string dbg_path = "";
    std::string fasta_path = "";
    std::string bam_path = "";
    std::string region = "";
    std::vector<std::string> subtract_v;
    std::vector<std::string> remove_v;
    std::string tx_kmer_path = "";
    bool low_coverage_filter = false;
    std::string ref_filter_path = "";
    bool keep_canary_kmers = false;
    std::string canary_path = "";
    std::string od = "";
    std::unordered_map<std::string, std::string> bam_files;

    float noise_size = 0;

    // Argument parsing
    CLI::App app{"Generates a de Bruijn Graph from fastx/bam files and counts the abundances of the constituent unitigs."};
    app.add_flag("-v,--verbose", verbose, "Print information about the pruning process to stdout.");
    app.add_flag("-d,--dry_run", dry_run, "Constructs and prunes a DBG, and counts abundances without saving to disk.");
    app.add_option("-t,--threads", threads, "Maximum number of threads to be used.");
    app.add_option("-g,--graph", dbg_path, "An existing de Bruijn graph to be pruned.");
    app.add_option("-f,--fasta", fasta_path, "A single FASTX file containing RNA reads for association.");
    app.add_option("-b,--bam", bam_path, "A BAM file or list of BAM files containing RNA reads.");
    app.add_option("-n,--region", region, "Only used with --bam/-b, usage: --region chrN:lb-ub.");
    app.add_option("-s,--subtract", subtract_v, "DBGs to be subtracted from the target DBG.");
    app.add_option("-r,--remove-kmers", remove_v, "Newline separated kmers whose encapsulating unitigs are to be removed from the target DBG.");
    app.add_option("-l,--transcriptomic-filter", tx_kmer_path, "Assume unitigs that have coverage < max(1, 0.5\% of median individual transcriptomic coverage) for an individual are sequencing errors and prune them. arg is a newline separated list of the kmers in all annotated transcripts of the gene.");
    app.add_flag("-c,--low-coverage-filter", low_coverage_filter, "Filter out unitigs with per-individual coverage low w.r.t. the per-individual coverage of the H1 neighborhood around the unitig.");
    app.add_option("--ref-fasta-filter", ref_filter_path, "Fasta file containing one sequence per entry. We find all kmers within H1 of the constituent kmers and prune the low coverage ones.");
    app.add_option("--noise", noise_size, "Add/subtract Poisson-distributed noise to count");

    app.add_option("--canary", canary_path, "List of sanity check kmers to check how they are retained after the individual pruning steps.");
    app.add_flag("--keep", keep_canary_kmers, "Keep all counts for canary unitigs.");
    app.add_option("-o,--output", od, "Directory to which the output files are to be saved.");
    CLI11_PARSE(app, argc, argv);

    if (fasta_path == "" && bam_path == "") {
        std::cerr << "You must supply reads, either in a fasta file using the -f flag or in a list or BAM files using the -b flag" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (!dry_run && od == "") {
        std::cerr << "You must supply an output directory using the -o flag" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Print the command used to run this program to od/run_info
    if (!dry_run) {
        ofstream run_info(od + "/run_info");
        for (int i = 0; i < argc; ++i) {
            run_info << argv[i] << " ";
        }
        run_info.close();
    }

    // Build graph
    ACDBG_Build_opt opt;
    opt.verbose = verbose;
    opt.debug = false;
    opt.nb_threads = threads;

    AbundanceCDBG t_graph;
    if (dbg_path != "") {
        std::cout << "Reading graph from file: " << dbg_path << std::endl;
        t_graph.read(dbg_path, opt.verbose);
    } else if (bam_path != "") {
        // Construct graph from BAM file or list of BAM files and prune it
        if (region == "") {
            std::cerr << "Parameter --region must be supplied when constructing graph from BAM files." << std::endl;
            exit(1);
        }
        if (opt.verbose) std::cout << "Constructing graph from BAM files: " << bam_path << std::endl;
        std::string id, path;

        // Check whether a single BAM file or list of BAM files
        std::ifstream infile(bam_path);
        infile >> id >> path;

        const size_t last_point = path.find_last_of(".");
        std::string path_ext = path.substr(last_point + 1);
        if (path_ext == "bam") {
            bam_files[id.substr(0, 7)] = path;
            while (infile >> id >> path) {
                bam_files[id.substr(0, 7)] = path;
            }
        } else {
            bam_files["PNXXXXX"] = bam_path;
        }

        size_t n_cohort = bam_files.size(), idx = 0;
        std::vector<std::string> bam_paths;
        bam_paths.reserve(n_cohort);
        t_graph.pns.reserve(n_cohort);
        t_graph.pn2idx.reserve(n_cohort);
        for (const auto &elem : bam_files) {
            t_graph.pns.push_back(elem.first);
            t_graph.pn2idx[elem.first] = idx;
            bam_paths.push_back(elem.second);
            ++idx;
        }

        /*
        std::for_each(bam_files.begin(), bam_files.end(), [&](std::pair<const std::string, std::string>  &elem){
            bam_paths.push_back(elem.second);
        });
        */
        if (opt.verbose) std::cout << "PruneACDBG: Constructing graph from " << bam_paths.size() << " BAM files" << std::endl;

        t_graph.build(opt, bam_paths, region, opt.verbose);

    } else {
        // Construct graph from fasta file and prune it

        // Pass through FASTX file in order to get pn2idx and idx2pn maps
        std::string line, pn, new_pn;
        size_t idx;
        std::ifstream infile(fasta_path);
        infile >> line;
        pn = line.substr(1, LEN_PN);
        idx = t_graph.pns.size();
        t_graph.pns.push_back(pn);
        t_graph.pn2idx[pn] = idx;
        while (infile >> line) {
            if (line.substr(0, 1) == ">") {
                new_pn = line.substr(1, LEN_PN);
                if (pn != new_pn){
                    pn = new_pn;
                    t_graph.pns.push_back(pn);
                    t_graph.pn2idx[pn] = ++idx;
                }
            }
        }

        if (opt.verbose) std::cout << "PruneACDBG: Constructing graph from fasta file: " << fasta_path << std::endl;
        t_graph.build(opt, fasta_path, opt.verbose);

    }

    std::vector<Kmer> canary;
    if (canary_path != "") {
        std::ifstream cinfile(canary_path);
        std::string seq;
        while (cinfile >> seq) {
            canary.push_back(Kmer(seq.c_str()));
        }
    }


    if (opt.debug) t_graph.write(od + "/unpruned_graph", opt.nb_threads, opt.verbose);

    // t_graph.flag_remove_tips(opt.verbose);

    // ==============
    // Pruning driver
    // ==============

    // Graph subtractions
    for (const auto &graph_path : subtract_v) {
        t_graph.flag_subtract_graph(graph_path, opt.nb_threads, opt.verbose);
    }

    // Kmer subtractions
    for (const auto &kmer_path : remove_v) {
        t_graph.flag_subtract_kmers(kmer_path, opt.verbose);
    }

    t_graph.reload_unflagged(opt, od, opt.verbose);
    if (canary.size() > 0) t_graph.canary_kmer_retention(canary);

    // Resize the abundance vectors in all the UnitigMaps in the graph to be
    // the size of the cohort.
    if (opt.verbose) std::cout << "PruneACDBG: Pre-count pruning done" << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    // Count kmers and attach to nodes
    if (bam_path == "") t_graph.count_kmers(fasta_path, LEN_PN, opt.verbose);
    else t_graph.count_kmers(bam_files, region, opt.verbose);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    if (opt.verbose) std::cout << "PruneACDBG: Kmer count took " << duration << " microseconds." << std::endl;

    t_graph.avg_counts(canary, keep_canary_kmers, true, opt.verbose);
    if (canary.size() > 0) t_graph.canary_kmer_retention(canary);

    if (noise_size > 0) {
        if (opt.verbose) std::cout << "Adding Poisson-distributed noise with parameter mu*" << noise_size << " to counts." << std::endl;
        t_graph.add_noise(noise_size);
    }

    // Prune w.r.t. individual transcriptomic median abundance.
    if (tx_kmer_path != "") {
        if (opt.verbose) std::cout << "Pruning w.r.t. individual transcriptomic median abundance." << std::endl;
        t_graph.prune_individual_mean(tx_kmer_path, true, canary, keep_canary_kmers, opt.verbose);
    }
    if (canary.size() > 0) t_graph.canary_kmer_retention(canary);

    // Prune low coverage kmers with high abundance kmers within H1 neighborhood.
    if (low_coverage_filter) {
        if (opt.verbose) std::cout << "Pruning w.r.t. abundance of unitigs in H1 neighborhood." << std::endl;
        t_graph.remove_low_wrt_h1(t_graph.pns.size()*0.05, opt.ratio_threshold, true, canary, keep_canary_kmers, opt.verbose);
    }
    if (canary.size() > 0) t_graph.canary_kmer_retention(canary);

    // Remove tips with cardinality > 1
    t_graph.flag_remove_tips(opt.verbose);
    if (canary.size() > 0) t_graph.canary_kmer_retention(canary);
    t_graph.remove_flagged_unitigs(canary, keep_canary_kmers, opt.verbose);
    if (canary.size() > 0) t_graph.canary_kmer_retention(canary);

    if (canary.size() > 0) t_graph.canary_kmer_retention(canary);

    // Calculate total number of reads per individual
    std::vector<float> idx2total;
    t_graph.tally_counts(idx2total);
    if (canary.size() > 0) t_graph.canary_kmer_retention(canary);

    if (!dry_run) {
        t_graph.write(od + "/graph", opt.nb_threads, opt.verbose);
        t_graph.write_abundances(od + "/pn2count");
        t_graph.write_abundances(od + "/pn2count.normalized", idx2total);
        t_graph.write_kmer2unitig(od + "/kmer2unitig");
    }
    return 0;
}
