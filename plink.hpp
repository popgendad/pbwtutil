#pragma once

#include <string.h>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <map>

namespace plink
{
// Unknown allele value in .bim file
const char UNKNOWN_VARIANT = '0';

// Two header bytes for .bed data
const int BED_MAGIC1 = 108;
const int BED_MAGIC2 = 27;

// Allele encodings
const unsigned char BED_HOMOZYGOUS_0 = (unsigned char)0;
const unsigned char BED_MISSING = (unsigned char)1;      // (see above)
const unsigned char BED_HETEROZYGOUS = (unsigned char)2; // (see above)
const unsigned char BED_HOMOZYGOUS_1 = (unsigned char)3;

// or ".hap" data if we have a phased interpretation
const int HAP_MAGIC1 = BED_MAGIC1;
const int HAP_MAGIC2 = 211;
const unsigned char HAP_00 = (unsigned char)0;
const unsigned char HAP_01 = (unsigned char)1;
const unsigned char HAP_10 = (unsigned char)2;
const unsigned char HAP_11 = (unsigned char)3;

// One header byte for matrix data orientation
const int SNP_MAJOR_ORDER = 1;        // SNP-major order
const int INDIVIDUAL_MAJOR_ORDER = 0; // individual-major order

// Define a snp_t type
typedef struct
{
    int ch; // chromosome number--kept as int for sorting purposes
    std::string rs;   // marker ID--generally matches '^rs[0-9]+$'
    double cM;        // recombination distance from previous marker (or beginning of chromosome)
    unsigned long bp; // coordinate location in base pairs
    char a0, a1;      // e.g., 'A','G'; represented by 0 and 1 in a .bed (or .hap) file, respectively (typically allele0 would be the major allele, but this may not be guaranteed).  Set to ALLELE_MISSING if the allele value is as of yet unknown.
} snp_t;

// Array of SNP data
typedef std::vector<snp_t> bim_t;

// Read one SNP from file, tell if you expect bim (6 column) or not (4 column, map) input
snp_t load_snp(std::istream &in);

// Read file into array
bim_t load_bim(std::istream &in);

// Write array to file
void dump_bim(std::ostream &, const bim_t &bim);

// Individual entry from FAM file
typedef struct
{
    std::string fid, iid, pid, mid, sex, phe;
} indiv_t;
typedef std::vector<indiv_t> fam_t;

// Individual entry from REG file
typedef struct
{
    std::string fid, iid, pid, mid, sex, phe, pop, reg;
} refsam_t;

typedef std::vector<refsam_t> reg_t;

indiv_t load_indiv(std::istream &in);

refsam_t load_refsam(std::istream &in);

void dump_indiv(std::ostream &out, const indiv_t &indiv);

void dump_refsam(std::ostream &out, const refsam_t &refsam);

void dump_fam(std::ostream &, const fam_t &);

void dump_reg(std::ostream &, const reg_t &);

fam_t load_fam(std::istream &);

reg_t load_reg(std::istream &);

// Compute bed data binary record size based on number of elements in each record
// (e.g., if using SNP-major order, this is the number of inidiviuals)
inline size_t bed_record_size(int n_minor) { return n_minor / 4 + ((n_minor % 4) ? 1 : 0); }

// retreive 2-bit encoding for the given major (e.g., marker index if we're doing SNP-major order) and minor (e.g., individual index) matrix entry
inline unsigned char get_encoding(const unsigned char *data, size_t record_size, size_t major, size_t minor)
{
    return (data[major * record_size + minor / 4] >> 2 * (minor % 4)) & 3; // mask = 3 = 00000011
}

// Update the encoding with a new two-bit value (e.g., BED_HOMOZYGOUS_0)
void set_encoding(unsigned char *data, size_t record_size, size_t major, size_t minor, const unsigned char genotype);

// Binary genotype data and associated functionality
class bed_t
{

public:
    const int header1, header2; // header bytes
    const int orientation;
    const size_t record_size; // number of bytes per record (major-order singleton array)
    const size_t size;        // number of bytes of data
    unsigned char *data;

    // create bed
    bed_t(int h1, int h2, int snp_major_order, size_t rsz, size_t dsz, unsigned char *a) : header1(h1), header2(h2), orientation(snp_major_order), record_size(rsz), size(dsz), data(a) {}
    bed_t(int snp_major_order, size_t rsz, size_t dsz, unsigned char *a) : header1(BED_MAGIC1), header2(BED_MAGIC2), orientation(snp_major_order), record_size(rsz), size(dsz), data(a) {}

    inline bool phased() const { return header2 == HAP_MAGIC2; }
    inline unsigned char genotype(size_t i, size_t m) { return get_encoding(data, record_size, orientation ? m : i, orientation ? i : m); }
    inline void set_genotype(size_t i, size_t m, unsigned char g) { set_encoding(data, record_size, orientation ? m : i, orientation ? i : m, g); }
};

// Create a new bed object and dynamically allocate data for it
bed_t *allocate_bed(size_t n_indiv, size_t n_snps, int orientation, bool phased = false);

// Load bed.  Use given unsigned char array or allocate a new one if argument is NULL
bed_t *load_bed(FILE *fin, size_t n_indiv, size_t n_snps, unsigned char *data = NULL);

// Dump bed, return number of bytes written
size_t dump_bed(FILE *, const bed_t *);

// Main plink_t class definition
class plink_t
{
public:
    bed_t *bed;
    bim_t bim;
    fam_t fam;

    // Access
    inline unsigned char genotype(size_t i, size_t m) const { return bed->genotype(i, m); }

    // If you know the indiv. and SNP but don't care about the orientation
    inline void set_genotype(size_t i, size_t m, const unsigned char g) { bed->set_genotype(i, m, g); }

    // For phased interpretation:
    inline bool haplotype(size_t i, size_t m, int parent) const { return (this->genotype(i, m) & (1 << parent)) ? true : false; }
    void set_haplotype(size_t indiv_index, size_t snp_index, int parent, bool bit);

    // Create custom
    plink_t() : bed(NULL) {}

};

std::ostream &operator<<(std::ostream &out, const plink_t &dataset);

// exception for I/O methods
class plink_exception_t : public std::exception
{

public:
    std::ostringstream oss;
    plink_exception_t() {}
    plink_exception_t(const plink_exception_t &peer) { oss << peer.what(); }
    plink_exception_t(const std::string &msg) { oss << msg; }
    template <typename P>
    plink_exception_t &operator<<(const P &msg)
    {
        oss << msg;
        oss.flush();
        return *this;
    }
    const char *what() const throw()
    {
        char *msg = new char[oss.str().size() + 1];
        strcpy(msg, oss.str().c_str());
        return msg; // not sure why returning oss.str().c_str() isn't good enough, but the error message is not always shown
        // (leak)
    }
    ~plink_exception_t() throw() {}
};

// Index a BIM
std::map<std::string, size_t> bimindex(const bim_t &bim);

// Index a FAM
std::map<std::string, size_t> famindex(const fam_t &fam);

// Index a REG
std::map<std::string, std::string> regindex(const reg_t &reg);

// Load available data
plink_t load(const std::string &stub);

// Dump data to disk
void dump(const plink_t &dataset, const std::string &bedfile, const std::string &bimfile, const std::string &famfile);

} //namespace
