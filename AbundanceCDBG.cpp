#include "AbundanceCDBG.hpp"

const int K = 31;

AbundanceCDBG::AbundanceCDBG() {}

bool AbundanceCDBG::read(const std::string &dbg_path, bool verbose=false) {
    return dbg.read(dbg_path, verbose);
}

bool AbundanceCDBG::build(const ACDBG_Build_opt &opt,
                          const std::vector<std::string> &bam_paths,
                          const std::string &region,
                          bool verbose=false) {
    return false;
    /*
    CDBG_Build_opt c_opt;
    c_opt.k = K;
    c_opt.nb_threads = opt.nb_threads;
    c_opt.build = true;
    c_opt.clipTips = true;
    c_opt.deleteIsolated = true;
    c_opt.region = region;
    c_opt.verbose = opt.verbose;
    c_opt.filename_seq_in = bam_paths;

    bool build_finished = dbg.build(c_opt);
    return build_finished && dbg.simplify(c_opt.deleteIsolated, c_opt.clipTips, c_opt.verbose);
    */
}

bool AbundanceCDBG::build(const ACDBG_Build_opt &opt,
                          const std::string &fasta_path,
                          bool verbose=false) {
    CDBG_Build_opt c_opt;
    c_opt.k = K;
    c_opt.nb_threads = opt.nb_threads;
    c_opt.build = true;
    // c_opt.clipTips = true;
    // c_opt.deleteIsolated = true;
    c_opt.clipTips = false;
    c_opt.deleteIsolated = false;
    c_opt.verbose = opt.verbose;
    c_opt.filename_seq_in.push_back(fasta_path);

    bool build_finished = dbg.build(c_opt);
    return build_finished && dbg.simplify(c_opt.deleteIsolated, c_opt.clipTips, c_opt.verbose);
}

bool AbundanceCDBG::write(const std::string path, size_t nb_threads=1, bool verbose=false) {
    return dbg.write(path, nb_threads, true, verbose);
}

bool AbundanceCDBG::remove_kmers(std::vector<std::string> &kmers_to_remove,
                                 bool verbose=false) {
    bool success = 1;
    Kmer kmer;
    UnitigMap<Node> um;
    for (auto const &seq : kmers_to_remove) {
        kmer = Kmer(seq.c_str());
        um = dbg.find(kmer);
        if (!um.isEmpty) {
            success = success && dbg.remove(um);
        }
    }
    if (verbose) std::cout << "AbundanceCDBG::remove_kmers(): Removed " << kmers_to_remove.size() << " unitigs." << std::endl;
    if (verbose) std::cout << "AbundanceCDBG::remove_kmers(): After " << dbg.size() << " unitigs." << std::endl;
    return success;
}

bool AbundanceCDBG::remove_kmers(std::vector<Kmer> &kmers_to_remove, bool verbose=false) {
    bool success = 1;
    UnitigMap<Node> um;
    for (auto const &kmer : kmers_to_remove) {
        um = dbg.find(kmer);
        if (!um.isEmpty) {
            success = success && dbg.remove(um);
        }
    }
    if (verbose) std::cout << "AbundanceCDBG::remove_kmers(): Removed " << kmers_to_remove.size() << " unitigs." << std::endl;
    if (verbose) std::cout << "AbundanceCDBG::remove_kmers(): After " << dbg.size() << " unitigs." << std::endl;
    return success;
}


void AbundanceCDBG::flag_kmers(std::vector<Kmer> &kmers_to_flag, bool verbose=false) {
    UnitigMap<Node> um;
    int flagged = 0;
    for (const auto &kmer : kmers_to_flag) {
        um = dbg.find(kmer);
        if (!um.isEmpty) {
            um.getData()->flag(um.dist, um.dist+um.len);
            ++flagged;
        }
    }
    if (verbose) {
        std::cout << "AbundanceCDBG::flag_kmers(): Flagged " << flagged << " kmers." << std::endl;
    }
}

void AbundanceCDBG::flag_unitigs(std::vector<Kmer> &unitigs_to_flag, bool verbose=false) {
    UnitigMap<Node> um;
    int flagged = 0;
    for (const auto &kmer : unitigs_to_flag) {
        um = dbg.find(kmer);
        if (!um.isEmpty) {
            um.getData()->flag(0, um.size-30);
            ++flagged;
        }
    }
    if (verbose) {
        std::cout << "AbundanceCDBG::flag_unitigs(): Flagged " << flagged << " unitigs." << std::endl;
    }
}

bool AbundanceCDBG::remove_flagged_unitig(UnitigMap<Node> &um,
                                          std::vector<Kmer> &canary,
                                          bool keep, bool verbose=false) {

    // BEGIN DEBUG
    verbose = false;
    // END DEBUG

    Node *data = um.getData();
    if (um.isEmpty || data->r.isEmpty()) {
        // Unitig is either nonextant or extant and not flagged; keep entire unitig
        return true;
    } else if (data->r.cardinality() == um.size-30) {
        // Entire unitig is tagged; delete entire unitig
        return dbg.remove(um, verbose);
    }

    std::vector<std::pair<uint32_t, float> > elems;
    std::vector<std::pair<uint32_t, float> > new_elems;
    SparseVector<float> abundance(data->abundance);
    abundance.getElements(elems);
    float read_counts = data->read_counts;

    std::string seq = um.referenceUnitigToString();
    // BEGIN DEBUG
    if (is_in_kmer_list(um, canary) || um.getUnitigHead().toString() == "ACGGCTGCCCGAAGCCCCCCGAGATTGCACT") {
        std::cout << "============================" << std::endl;
        std::cout << "Deleting an important unitig" << std::endl;
        std::cout << "============================" << std::endl;
    }
    // END DEBUG
    Kmer kmer;
    UnitigMap<Node> new_um;
    // Get array of the indices of contiguous unflagged kmers from bit
    // vector, i.e. the kmers we want to keep
    Roaring r = Roaring(data->r);

    r.flip(0, um.size-30);
    // We've stored all the data corresponding to um, so we can safely remove it
    dbg.remove(um, verbose);
    size_t n_flags = r.cardinality();
    uint32_t *arr = new uint32_t[n_flags];
    r.toUint32Array(arr);
    uint32_t lb = arr[0], ub = arr[0], len;
    bool ret = true;
    int created = 0;

    // Assert that we're working on the forward strand
    if (!um.strand) {
        // Reverse r[0, um.size-30]
        Roaring nr;
        for (int i = 0; i < um.size-30; ++i) {
            if (r.contains(i)) {
                nr.add(um.size-31-i);
            }
        }
        r = nr;
    }

    // DEBUG
    bool is_jct = 0;
    // END DEBUG

    // Create new unitigs from contiguous blocks in old unitig
    for (size_t i = 1; i < n_flags; ++i) {
        if (arr[i] == arr[i-1]+1) {
            ub = arr[i];
        } else {
            // Add new unitig from contiguous block
            len = ub - lb + 1;
            ret = ret && dbg.add(seq.substr(lb, len+30), verbose);
            ++created;
            kmer = Kmer(seq.substr(lb, 31).c_str());
            new_um = dbg.find(kmer);
            float w_new = len / (new_um.size - 30); // Unitig size is in characters, not kmers
            data = new_um.getData();
            data->abundance.getElements(new_elems);

            // DEBUG
            float before = 0., after = 0.;
            // END DEBUG

            // If the unitig was joined to an extant unitig, the abundances are
            // the weighted averages of the abundances of the new and old unitigs
            for (const auto &it : new_elems) {
                // DEBUG
                before += data->abundance[it.first];
                // END DEBUG
                data->abundance[it.first] = it.second*(1-w_new);
            }
            for (const auto &it : elems) {
                if (data->abundance.contains(it.first)) {
                    data->abundance[it.first] += it.second*w_new;
                } else {
                    data->abundance.insert(it.first, it.second*w_new);
                }
                // DEBUG
                after += data->abundance[it.first];
                // END DEBUG
            }
            new_elems.clear();

            data->read_counts += read_counts*w_new;
            lb = arr[i];
            ub = lb;
        }
    }
    delete [] arr;

    // Last new unitig
    len = ub - lb + 1;
    ret = ret && dbg.add(seq.substr(lb, len+30), verbose);
    ++created;
    kmer = Kmer(seq.substr(lb, 31).c_str());
    new_um = dbg.find(kmer);
    float w_new = float(len) / (new_um.size - 30); // Unitig size is in characters, not kmers
    data = new_um.getData();
    data->abundance.getElements(new_elems);
    // If the unitig was joined to an extant unitig, the abundances are
    // the weighted averages of the abundances of the new and old unitigs
    for (const auto &it : new_elems) {
        data->abundance[it.first] = it.second*(1-w_new);
    }
    for (const auto &it : elems) {
        if (data->abundance.contains(it.first)) {
            data->abundance[it.first] += it.second*w_new;
        } else {
            data->abundance.insert(it.first, it.second*w_new);
        }
    }
    data->read_counts += read_counts*w_new;

    if (verbose) {
        std::cout << "AbundanceCDBG::remove_flagged_unitig(): ";
        std::cout << "Removed 1 unitig, created " << created << " unitigs." << std::endl;
    }

    return ret;
}

bool AbundanceCDBG::remove_flagged_unitigs(std::vector<Kmer> &canary, bool keep, bool verbose=false) {
    std::vector<Kmer> flagged_unitigs;

    for (const auto &um : dbg) {
        if (!um.isEmpty && !um.getData()->r.isEmpty()) {
            flagged_unitigs.push_back(um.getUnitigHead());
        }
    }
    if (verbose) std::cout << "AbundanceCDBG::remove_flagged_unitigs(): Removing " << flagged_unitigs.size() << " unitigs." << std::endl;

    UnitigMap<Node> um2rm;
    for (auto &kmer : flagged_unitigs) {
        um2rm = dbg.find(kmer);
        remove_flagged_unitig(um2rm, canary, keep, verbose);
    }
    return true;
}

bool AbundanceCDBG::subtract_graph(const std::string &path, size_t threads=1, bool verbose=false) {
    CompactedDBG<Node> m_dbg;
    if (verbose) std::cout << "AbundanceCDBG::subtract_graph(): Reading graph from file " << path << std::endl;
    m_dbg.read(path, threads);
    bool success = 1;
    int out = 0;

    // BEGIN hack to bypass Bifrost bug
    std::string troublesome_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG";
    Kmer troublesome_kmer(troublesome_seq.c_str());
    UnitigMap<Node> troublesome_unitig = dbg.find(troublesome_kmer);
    if (!troublesome_unitig.isEmpty) {
        dbg.remove(troublesome_unitig);
    }
    // END hack to bypass Bifrost bug

    for (auto& m_um : m_dbg) {
        std::string unitig = m_um.mappedSequenceToString();
        KmerIterator it_end, it(unitig.c_str());
        for (; it != it_end; ++it) {
            UnitigMap<Node> um = dbg.find(it->first);
            if (!um.isEmpty) {
                success = success && dbg.remove(um);
                if (success) ++out;
            }
        }
    }

    if (verbose) {
        if (success) {
            std::cout << "AbundanceCDBG::subtract_graph(): ";
            std::cout << "Removed " << out << " unitigs." << std::endl;
            std::cout << "AbundanceCDBG::subtract_graph(): ";
            std::cout << "After: " << dbg.size() << " unitigs." << std::endl;
        } else {
            std::cout << "Error removing one or more unitigs." << std::endl;
        }
    }
    return success;
}

bool AbundanceCDBG::subtract_kmers(const std::string &path, bool verbose=false) {
    ifstream infile(path);
    std::string seq;
    Kmer kmer;
    UnitigMap<Node> um;
    bool success = 1;
    int out = 0;
    while (infile >> seq) {
        kmer = Kmer(seq.c_str());
        um = dbg.find(kmer);
        if (!um.isEmpty) {
            success = success && dbg.remove(um);
            if (success) ++out;
        }
    }
    if (verbose) {
        if (success) {
            std::cout << "AbundanceCDBG::subtract_kmers(): ";
            std::cout << "Removed " << out << " unitigs." << std::endl;
            std::cout << "AbundanceCDBG::subtract_kmers(): ";
            std::cout << "After: " << dbg.size() << " unitigs." << std::endl;
        } else {
            std::cout << "Error removing one or more unitig" << std::endl;
        }
    }
    return success;
}

bool AbundanceCDBG::flag_subtract_graph(const std::string &path, size_t threads=1, bool verbose=false) {
    CompactedDBG<Node> m_dbg;
    std::cout << "AbundanceCDBG::subtract_graph(): Reading graph from file " << path << std::endl;
    m_dbg.read(path, threads);
    int flagged = 0;

    for (auto& m_um : m_dbg) {
        std::string unitig = m_um.mappedSequenceToString();
        KmerIterator it_end, it(unitig.c_str());
        for (; it != it_end; ++it) {
            UnitigMap<Node> um = dbg.find(it->first);
            if (!um.isEmpty) {
                um.getData()->flag(um);
                ++flagged;
            }
        }
    }
    if (verbose) {
        std::cout << "AbundanceCDBG::subtract_graph(): ";
        std::cout << "Flagged " << flagged << " unitigs for removal." << std::endl;
    }
    return true;
}

bool AbundanceCDBG::flag_subtract_kmers(const std::string &path, bool verbose=false) {
    ifstream infile(path);
    std::string seq;
    Kmer kmer;
    UnitigMap<Node> um;
    int flagged = 0;
    while (infile >> seq) {
        kmer = Kmer(seq.c_str());
        um = dbg.find(kmer);
        if (!um.isEmpty) {
            um.getData()->flag(um);
            ++flagged;
        }
    }
    if (verbose) {
        std::cout << "AbundanceCDBG::subtract_kmers(): ";
        std::cout << "Flagged " << flagged << " unitigs for removal." << std::endl;
    }
    return true;
}

// Writes sequences of unflagged unitigs to a temporary fasta file, reads them
// into a new graph and deletes the temporary file.
bool AbundanceCDBG::reload_unflagged(ACDBG_Build_opt &opt, std::string od="", bool verbose=false) {
    std::string filename = std::string(std::tmpnam(nullptr));
    size_t slash = filename.find_last_of("/");
    filename = od + filename.substr(slash) + ".fa";
    size_t id = 0;
    if (verbose) std::cout << "AbundanceCDBG::reload_unflagged(): Writing sequence to temporary file " << filename << std::endl;
    ofstream outfile(filename);
    for (const auto &um : dbg) {
        if (um.getData()->r.isEmpty()) {
            outfile << ">" + std::to_string(id++) + "\n" << um.referenceUnitigToString() << "\n";
        }
    }
    outfile.close();

    CDBG_Build_opt c_opt;
    c_opt.k = K;
    c_opt.nb_threads = opt.nb_threads;
    c_opt.build = true;
    c_opt.clipTips = true;
    c_opt.deleteIsolated = true;
    c_opt.verbose = opt.verbose;
    c_opt.filename_ref_in.push_back(filename);

    // CompactedDBG::clear() sets the invalid flag to true, so we have to create
    // a new object
    dbg = CompactedDBG<Node>();
    bool ret = dbg.build(c_opt);
    ret = ret && dbg.simplify(c_opt.deleteIsolated, c_opt.clipTips, c_opt.verbose);

    if (verbose) std::cout << "AbundanceCDBG::reload_unflagged(): Deleting temporary file " << filename << std::endl;
    std::remove(filename.c_str());
    if (verbose) std::cout << "AbundanceCDBG::reload_unflagged(): Deleted temporary file " << filename << std::endl;

    return ret;
}

// Not normalized
void AbundanceCDBG::write_abundances(std::string path) {
    ofstream outfile(path);
    Node* data;
    Kmer kmer;
    std::vector<std::pair<uint32_t, float> > elems;
    for (const auto &um : dbg) {
        data = um.getData();
        kmer = um.getUnitigHead();
        data->abundance.getElements(elems);
        for (const auto &it : elems) {
            if (it.second > 0) outfile << pns[it.first] << "\t" << kmer.toString() << "\t" << it.second << std::endl;
        }
        elems.clear();
    }
    outfile.close();
}

// Normalized
void AbundanceCDBG::write_abundances(std::string path, std::vector<float> &idx2total) {
    ofstream outfile(path);
    Node* data;
    Kmer kmer;
    std::vector<std::pair<uint32_t, float> > elems;
    for (const auto &um : dbg) {
        data = um.getData();
        kmer = um.getUnitigHead();
        data->abundance.getElements(elems);
        for (const auto &it : elems) {
            if (it.second > 0) outfile << pns[it.first] << "\t" << kmer.toString() << "\t" << it.second/idx2total[it.first] << std::endl;
        }
        elems.clear();
    }
    outfile.close();
}

void AbundanceCDBG::write_kmer2unitig(std::string path) {
    ofstream outfile(path);
    for (const auto& um : dbg) {
        outfile << um.getUnitigHead().toString() << "\t" << um.referenceUnitigToString() << std::endl;
    }
}

void AbundanceCDBG::count_kmers(const std::string &fasta_path, size_t len_pn=7, bool verbose=false) {
    if (verbose) std::cout << "AbundanceCDBG::count_kmers(): Counting kmers from " << fasta_path << std::endl;
    std::ifstream infile(fasta_path);
    std::string line, pn, new_pn;
    size_t idx;
    UnitigMap<Node> um;
    while (std::getline(infile, line)) {
        if (line.substr(0, 1) == ">") {
            new_pn = line.substr(1, len_pn);
            if (new_pn != pn) {
                pn = new_pn;
                idx = pn2idx[pn];
            }
        } else {
            // Skip reads that are shorter than K (if any)
            if (line.length() < K) continue;
            int i = 0; // Number of kmers from line found
            while (i < line.length() - K + 1) {
                um = dbg.findUnitig(line.c_str(), i, line.length());
                if (!um.isEmpty) {
                    i += um.len;
                    Node* n = um.getData();
                    n->read_counts += um.len;
                    // Add one count for each kmer in the mapping
                    if (n->abundance.contains(idx)) {
                        n->abundance[idx] += um.len;
                    } else {
                        n->abundance.insert(idx, um.len);
                    }
                } else {
                    ++i;
                }
            }
        }
    }
}

void AbundanceCDBG::count_kmers(const std::unordered_map<std::string, std::string> &bam_files,
                                const std::string &region, bool verbose=false) {
    std::cerr << "AbundanceCDBG::count_kmers: BAM files not supported" << std::endl;
    exit(1);

    /*
    if (verbose) std::cout << "AbundanceCDBG::count_kmers(): Counting kmers from BAM files" << std::endl;
    BAMReader br;
    std::string read, pn;
    UnitigMap<Node> um;
    for (auto &it : bam_files) {
        BAMReader br = BAMReader(it.second, region);
        size_t idx = pn2idx[it.first];
        br.next();
        while (!br.is_empty()) {
            read = br.decode();
            int i = 0; // Number of kmers from the read that have been found
            while (i < read.length() - K + 1) {
                um = dbg.findUnitig(read.c_str(), i, read.length());
                if (!um.isEmpty) {
                    i += um.len;
                    Node* n = um.getData();
                    n->read_counts += um.len;
                    // Add one count for each kmer in the UnitigMap
                    if (n->abundance.contains(idx)) {
                        n->abundance[idx] += um.len;
                    } else {
                        n->abundance.insert(idx, um.len);
                    }
                } else {
                    ++i;
                }
            }
            br.next();
        }
        br.close();
    }
    */
}

void AbundanceCDBG::avg_counts(std::vector<Kmer> &canary, bool keep, bool flag=false, bool verbose=false) {
    if (verbose) std::cout << "AbundanceCDBG::avg_counts(): Setting unitig counts to mean kmer count." << std::endl;
    Kmer kmer;
    // The abundance of a unitig should be the average of the abundances of the
    // constituent kmers
    std::vector<std::pair<uint32_t, float> > elems;
    std::vector<Kmer> kmers_to_remove;
    for (auto const &um : dbg) {
        Node *data = um.getData();
        size_t n_kmers = um.len;
        data->abundance.getElements(elems);
        for (const auto &it : elems) {
            if (!(keep && is_in_kmer_list(um, canary)) && it.second == 1) {
                data->abundance.remove(it.first);
                data->read_counts -= 1;
            } else {
                data->abundance[it.first] /= n_kmers;
            }
        }

        if (data->read_counts > 0) {
            data->read_counts /= n_kmers;
        } else {
            kmers_to_remove.push_back(um.getUnitigHead());
        }
        elems.clear();
    }
    if (!flag) remove_kmers(kmers_to_remove, verbose);
    else flag_unitigs(kmers_to_remove, verbose);
}

void AbundanceCDBG::detect_tips(std::vector<Kmer> &tips) {
    for (auto& um : dbg) {
        if (is_tip(um)) {
            tips.push_back(um.getUnitigHead());
        }
    }
}

void AbundanceCDBG::remove_tips(bool verbose=false) {
    std::vector<Kmer> tips;
    if (verbose) std::cout << "AbundanceCDBG::remove_tips(): Detecting tips." << std::endl;
    detect_tips(tips);
    if (verbose) std::cout << "AbundanceCDBG::remove_tips(): Removing " << tips.size() << " tips." << std::endl;
    remove_kmers(tips, verbose);
}

void AbundanceCDBG::flag_remove_tips(bool verbose=false) {
    std::vector<Kmer> tips;
    if (verbose) std::cout << "AbundanceCDBG::remove_tips(): Detecting tips." << std::endl;
    detect_tips(tips);
    if (verbose) std::cout << "AbundanceCDBG::remove_tips(): Removing " << tips.size() << " tips." << std::endl;
    flag_unitigs(tips, verbose);
}

void AbundanceCDBG::add_noise(float size) {
    // Calculate each PN's mean expression over all nodes

    Node* data;
    std::vector<std::pair<uint32_t, float> > elems;
    std::vector<float> idx2mean;
    std::vector<uint32_t> lens;
    idx2mean.resize(pns.size(), 0.);
    lens.resize(pns.size(), 0);
    for (const auto um : dbg) {
        data = um.getData();
        data->abundance.getElements(elems);
        for (const auto &it : elems) {
            if (it.second > 0) {
                idx2mean[it.first] += it.second;
                ++lens[it.first];
            }
        }
        elems.clear();
    }

    float mean = 0.;
    for (size_t i = 0; i < idx2mean.size(); ++i) {
        mean += idx2mean[i] / lens[i];
    }
    mean = mean / idx2mean.size();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::poisson_distribution<> d(mean * size);

    for (auto um : dbg) {
        data = um.getData();
        for (size_t i = 0; i < pns.size(); ++i) {
            // Sample how many reads to add/subtract
            uint32_t events = d(gen);
            if (events > 0) {
                // Insert key into sparse vector if it was not previously present
                if (!data->abundance.contains(i)) {
                    data->abundance.insert(i, 0);
                }
                // Flip coin to see if we add or remove reads
                if (rand() % 2 == 1) {
                    data->abundance[i] += events;
                } else {
                    float diff = data->abundance[i] - events;
                    data->abundance[i] = (diff > 0) ? diff : 0;
                }
            }
        }
    }
}

// path points to a file containing all kmers in the annotated transcriptome
// for the current region, one kmer per line.
bool AbundanceCDBG::prune_individual_mean(const std::string &path, bool flag, std::vector<Kmer> &canary, bool keep, bool verbose=false) {
    if (verbose) std::cout << "AbundanceCDBG::prune_individual_mean(): Pruning w.r.t. transcriptomic median expression." << std::endl;
    ifstream infile(path);
    UnitigMap<Node> um;
    Node* data;
    std::string seq;
    Kmer kmer;
    std::vector<std::pair<uint32_t, float> > elems;
    std::vector<std::vector<float> > counts;
    counts.resize(pns.size(), {});
    while (infile >> seq) {
        kmer = Kmer(seq.c_str());
        um = dbg.find(kmer);
        if (um.isEmpty) continue;
        data = um.getData();
        data->abundance.getElements(elems);
        for (const auto &it : elems) {
            if (it.second > 0) {
                counts[it.first].push_back(it.second);
            }
        }
        elems.clear();
    }
    infile.close();

    std::vector<float> idx2med;
    idx2med.reserve(counts.size());
    for (size_t i; i < counts.size(); ++i) {
        size_t l = counts[i].size();
        if (l == 0) {
            idx2med.push_back(0);
        } else if (l == 1) {
            idx2med.push_back(counts[i][0]);
        } else if (l > 1) {
            std::sort(counts[i].begin(), counts[i].end());
            float median = (l % 2 == 0) ? counts[i][l/2] : (float)(counts[i][l/2] + counts[i][l/2 + 1]) / 2.;
            idx2med.push_back(median);
        }
    }

    int out = 0;
    float threshold = .005;
    std::vector<Kmer> kmers_to_remove;
    for (auto& um : dbg) {

        if (keep && is_in_kmer_list(um, canary)) continue;

        Node* data = um.getData();
        data->abundance.getElements(elems);
        for (const auto &it : elems) {
            if (idx2med[it.first] == 0) continue;
            if (it.second < 2 || (it.second / idx2med[it.first]) < threshold) {
                data->read_counts -= it.second;
                data->abundance.remove(it.first);
            }
        }

        if (data->read_counts < .1) {
            kmers_to_remove.push_back(um.getUnitigHead());
            ++out;
        }
        elems.clear();
    }

    if (!flag) {
        if (verbose) std::cout << "AbundanceCDBG::prune_individual_mean(): Removing " << out << " unitigs" << std::endl;
        remove_kmers(kmers_to_remove, verbose);
        if (verbose) std::cout << "AbundanceCDBG::prune_individual_mean(): After " << dbg.size() << " unitigs" << std::endl;
    } else {
        if (verbose) std::cout << "AbundanceCDBG::prune_individual_mean(): Flagging " << out << " unitigs" << std::endl;
        flag_unitigs(kmers_to_remove, verbose);
    }

    return true;
}

// Removes unitigs with low expression relative to their H1 neighbors.
//bool AbundanceCDBG::remove_low_wrt_h1(size_t count_threshold=50, float ratio_threshold=.25, bool flag=false, std::vector<Kmer> &canary, bool keep, bool verbose=false) {
bool AbundanceCDBG::remove_low_wrt_h1(size_t count_threshold, float ratio_threshold, bool flag, std::vector<Kmer> &canary, bool keep, bool verbose=false) {
    Node *data, *data_alt;
    UnitigMap<Node> um_alt;
    std::vector<Kmer> kmers_to_remove;
    std::set<size_t> idx_to_remove;
    std::vector<std::pair<size_t, UnitigMap<Node> > > h1;
    std::vector<std::pair<uint32_t, float> > elems;

    int out = 0;
    // Find unitigs with low aggregate expression, iterate over the pns that
    // contribute to the coverage and find kmers that are in H1(unitig).
    // We then check to see if any of those unitigs is expressed by the same pn
    // and have a substantially higher coverage than the initial unitig. If so,
    // we remove the reads for that pn from the initial unitig.
    for (const auto &um : dbg) {

        if (keep && is_in_kmer_list(um, canary)) continue;

        if (!um.isEmpty) {
            data = um.getData();
            if (data->read_counts < count_threshold) {
                data->abundance.getElements(elems);
                std::string seq = um.mappedSequenceToString();
                h1 = dbg.searchSequence(seq, 0, 0, 0, 1, 0);

                for (const auto &it : h1) {
                    // Kmers are unique in the graph
                    if (!it.second.isEmpty) {
                        data_alt = it.second.getData();
                        for (const auto &elem_it : elems) {
                            if (data_alt->abundance.contains(elem_it.first) &&
                                ((elem_it.second / data_alt->abundance[elem_it.first]) < ratio_threshold)) {
                                idx_to_remove.insert(elem_it.first);
                            }
                        }
                    }
                }
                for (const auto &idx : idx_to_remove) {
                    if (data->abundance.contains(idx)) {
                        data->read_counts -= data->abundance[idx];
                        data->abundance.remove(idx);
                    }
                }
                if (data->read_counts <= 1.) {
                    Kmer k = um.getUnitigHead();
                    kmers_to_remove.push_back(k);
                    ++out;
                }
                idx_to_remove.clear();
                h1.clear();
                elems.clear();
            }
        }
    }

    if (!flag) {
        if (verbose) std::cout << "AbundanceCDBG::remove_low_wrt_h1(): Removing " << out << " unitigs" << std::endl;
        remove_kmers(kmers_to_remove, verbose);
        if (verbose) std::cout << "AbundanceCDBG::remove_low_wrt_h1(): After " << dbg.size() << " unitigs" << std::endl;
    } else {
        if (verbose) std::cout << "AbundanceCDBG::remove_low_wrt_h1(): Flagging " << out << " unitigs" << std::endl;
        flag_unitigs(kmers_to_remove, verbose);
    }

    return true;
}

bool AbundanceCDBG::simplify(bool verbose=false) {
    return dbg.simplify(1, 1, verbose);
}

void AbundanceCDBG::tally_counts(std::vector<float> &idx2total) {
    // Calculate total expression per individual in graph
    Node* data;
    std::vector<std::pair<uint32_t, float> > elems;
    // Initialize idx2total with zeros
    idx2total.resize(pns.size(), 0.);
    for (const auto &um : dbg) {
        data = um.getData();
        data->abundance.getElements(elems);
        for (const auto &it : elems) {
            idx2total[it.first] += it.second;
        }
        elems.clear();
    }
}

// Sanity check
float AbundanceCDBG::canary_kmer_retention(std::vector<Kmer> &kmers) {
    float ratio = 0.;
    float flagged = 0.;
    float count = 0.;
    std::vector<std::pair<uint32_t, float> > elems;
    UnitigMap<Node> um;
    // START DEBUG
    std::vector<std::string> canary_pns;
    for (const auto &kmer : kmers) {
        um = dbg.find(kmer);
        if (!um.isEmpty) {
            ratio += 1.;
            // Check for flagged kmers
            if (!um.getData()->r.isEmpty()) {
                flagged += 1.;
            }
            um.getData()->abundance.getElements(elems);
            for (const auto &it : elems) {
                count += it.second;
                canary_pns.push_back(pns[it.first]);
            }
            elems.clear();
        }
    }
    std::cout << "\033[1;40mAbundanceCDBG::count_present_canary_kmers(): ";
    std::cout << ratio << " of " << kmers.size() << " canary kmers present in graph.\033[0m" << std::endl;
    std::cout << "\033[1;40mAbundanceCDBG::count_present_canary_kmers(): ";
    std::cout << flagged << " canary kmers flagged for deletion.\033[0m" << std::endl;
    std::cout << "\033[1;40mAbundanceCDBG::count_present_canary_kmers(): ";
    std::cout << count/kmers.size() << " mean expression of canary kmers.\033[0m" << std::endl;

    return ratio / static_cast<float>(kmers.size());
}

// Rudimentary splict junction test
// Supply a vector of all kmers overlapping potential splice junctions
bool AbundanceCDBG::is_in_kmer_list(const UnitigMap<Node> &um, std::vector<Kmer> &kmers) {

    for (int i = 0; i < um.len; ++i) {
        for (const auto kmer : kmers) {
            if (um.getMappedKmer(i).toString() == kmer.toString()) {
                std::cout << "AbundanceCDBG::is_in_kmer_list(): Found splice junction" << std::endl;
                return 1;
            }
        }
    }
    return 0;
}
