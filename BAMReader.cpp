#include <iostream>

#include "BAMReader.hpp"

BAMReader::BAMReader() {}

BAMReader::BAMReader(std::string const& path, std::string const& region) {
    open(path, region);
}

BAMReader::BAMReader(std::string const& path, std::string const& _id, std::string const& region) {
    id = _id;
    open(path, region);
}

int BAMReader::open(std::string const& path, std::string const& region) {
    fp = htslib::hts_open(path.c_str(), "r");
    if (!fp) {
        std::cerr << "[BAMReader] ERROR: Couldn't read file " << path << std::endl;
        std::exit(1);
    }
    fp->bam_header = htslib::sam_hdr_read(fp);
    init_itr(path, region);
    return 0;
}

int BAMReader::init_itr(std::string const& path, std::string const& region) {
    fp->idx = htslib::sam_index_load(fp, path.c_str());
    itr = htslib::sam_itr_querys(fp->idx, fp->bam_header, region.c_str());
    last_rec = htslib::bam_init1();
    return (itr == nullptr);
}

int BAMReader::close() {
    htslib::hts_itr_destroy(itr);
    itr = nullptr;
    htslib::bam_destroy1(last_rec);
    last_rec = nullptr;
    htslib::hts_close(fp);
    fp = nullptr;
    return 0;
}

int BAMReader::next() {
    last_itr = htslib::sam_itr_next(fp, itr, last_rec);
}

std::string BAMReader::decode() const {
    std::string seq;
    htslib::bam1_core_t const & core = last_rec->core;
    seq.resize(core.l_qseq);
    uint8_t* seqptr bam_get_seq(last_rec);
    for (int32_t i = 0; i < core.l_qseq; ++i) {
        seq[i] = magic_string[bam_seqi(seqptr, i)];
    }
    return seq;
}

bool BAMReader::is_empty() const {
    return (last_itr < 0);
}

// Reworked logic from Hannes's GraphTyper
// lhs > rhs if the sequence starts later or starts at the same locus and is
// further back in the alphabet.
bool BAMReader::operator>(BAMReader const& rhs) const {
    if (last_rec->core.tid > rhs.last_rec->core.tid) return true;
    else if (rhs.last_rec->core.tid > last_rec->core.tid) return false;
    else if (last_rec->core.pos > rhs.last_rec->core.pos) return true;
    else if (rhs.last_rec->core.pos > last_rec->core.pos) return false;
    else if (last_rec->core.l_qseq > rhs.last_rec->core.l_qseq) return true;
    else if (rhs.last_rec->core.l_qseq > last_rec->core.l_qseq) return false;

    int32_t const l_qseq = (last_rec->core.l_qseq + 1) / 2;
    uint8_t* lhs_seqptr = bam_get_seq(last_rec);
    uint8_t* rhs_seqptr = bam_get_seq(rhs.last_rec);

    for (int32_t i = 0; i < l_qseq; ++i) {
        auto lhs_seq = *(lhs_seqptr + i);
        auto rhs_seq = *(rhs_seqptr + i);

        if (lhs_seq > rhs_seq) return true;
        else if (lhs_seq > rhs_seq) return false;
    }
    return false;
}

bool BAM_Parser::is_bam(std::string& path) {
    htslib::hFILE* hfile = htslib::hopen(path.c_str(), "r");
    htslib::htsFormat fmt;

    if (htslib::hts_detect_format(hfile, &fmt)) {
        std::cerr << "Could not open BAM file: " << path << std::endl;
        return false;
    }

    if (htslib::hclose(hfile)) {
        std::cerr << "Could not close BAM file: " << path << std::endl;
    }

    return fmt.format == htslib::bam;
}
