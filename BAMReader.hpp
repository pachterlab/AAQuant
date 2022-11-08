#ifndef BAMReader_h
#define BAMReader_h

#include <string>

namespace htslib {
    #include <htslib/hts.h>
    #include <htslib/hfile.h>
    #include <htslib/sam.h>
}

class BAMReader {
    private:
        htslib::htsFile* fp = nullptr;
        htslib::hts_itr_t* itr = nullptr;
        int last_itr = 0;
        std::string magic_string = "=ACMGRSVTWYHKDBN";

    public:
        std::string id = "";
        htslib::bam1_t* last_rec = nullptr;

        BAMReader();
        BAMReader(std::string const& path, std::string const& region);
        BAMReader(std::string const& path, std::string const& id, std::string const& region);
        int open(std::string const& path, std::string const& region);
        int close();
        int next();
        int init_itr(std::string const& path, std::string const& region);
        std::string decode() const;
        bool is_empty() const;
        bool operator>(BAMReader const& rhs) const;

        static bool is_bam(std::string& path);
};
#endif
