#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include <fstream>
#include <map>
#include <mutex>
#include "plink.hpp"

namespace plink
{

// mutex for doing io synchronized, used by load_bim and load_fam
std::mutex io_mutex;

// Read one SNP from input stream
snp_t load_snp(std::istream &in)
{
    snp_t snp;
    static std::string line;
    getline(in, line);
    std::istringstream iss(line);
    std::string ch;
    iss >> ch >> snp.rs >> snp.cM >> snp.bp >> std::ws;
    if (!iss.eof())
        iss >> snp.a0 >> std::ws;
    else
        snp.a0 = UNKNOWN_VARIANT;
    if (!iss.eof())
        iss >> snp.a1;
    else
        snp.a1 = UNKNOWN_VARIANT;

    // set snp.ch to an integer
    //     X    X chromosome                    -> 23
    //     Y    Y chromosome                    -> 24
    //     XY   Pseudo-autosomal region of X    -> 25
    //     MT   Mitochondrial                   -> 26
    if (ch == "X")
        snp.ch = 23;
    else if (ch == "Y")
        snp.ch = 24;
    else if (ch == "XY")
        snp.ch = 25;
    else if (ch == "MT")
        snp.ch = 26;
    else
    {
        std::istringstream ch_iss(ch);
        ch_iss >> snp.ch;
    }
    return snp;
}

// Load BIM file
bim_t load_bim(std::istream &in)
{
    bim_t bim;
    io_mutex.lock();
    for (in >> std::ws; !in.eof(); in >> std::ws)
        bim.push_back(load_snp(in));
    io_mutex.unlock();
    return bim;
}

// Write one SNP to output stream
void dump_snp(std::ostream &out, const snp_t &snp)
{
    static const char delim = '\t';
    out << snp.ch << delim << snp.rs << delim << snp.cM << delim << snp.bp;
    out << delim << snp.a0 << delim << snp.a1 << '\n';
}

// Write entire BIM file to output stream
void dump_bim(std::ostream &out, const bim_t &bim)
{
    for (size_t i = 0; i < bim.size(); ++i)
        dump_snp(out, bim[i]);
}

// Dump FAM sample to output stream
void dump_indiv(std::ostream &out, const indiv_t &indiv)
{
    static const char delim = '\t';
    out << indiv.fid << delim << indiv.iid << delim << indiv.pid << delim << indiv.mid << delim << indiv.sex << delim << indiv.phe;
}

// Dump REG sample to output stream
void dump_refsam(std::ostream &out, const refsam_t &refsam)
{
    static const char delim = '\t';
    out << refsam.fid << delim << refsam.iid << delim << refsam.pid << delim << refsam.mid << delim << refsam.sex;
    out << delim << refsam.phe << delim << refsam.pop << delim << refsam.reg;
}

// Dump FAM to output stream
void dump_fam(std::ostream &out, const fam_t &fam)
{
    for (size_t i = 0; i < fam.size(); ++i)
    {
        dump_indiv(out, fam[i]);
        out << '\n';
    }
}

// Dump REG to output stream
void dump_reg(std::ostream &out, const reg_t &reg)
{
    for (size_t i = 0; i < reg.size(); ++i)
    {
        dump_refsam(out, reg[i]);
        out << '\n';
    }
}

// Load individual from FAM
indiv_t load_indiv(std::istream &in)
{
    indiv_t indiv;
    in >> indiv.fid >> indiv.iid >> indiv.pid >> indiv.mid >> indiv.sex >> indiv.phe;
    return indiv;
}

// Load individual from REG
refsam_t load_refsam(std::istream &in)
{
    refsam_t refsam;
    in >> refsam.fid >> refsam.iid >> refsam.pid >> refsam.mid >> refsam.sex >> refsam.phe >> refsam.pop >> refsam.reg;
    return refsam;
}

// Load FAM file from stream
fam_t load_fam(std::istream &in)
{
    fam_t fam;
    io_mutex.lock();
    for (in >> std::ws; !in.eof(); in >> std::ws)
        fam.push_back(load_indiv(in));
    io_mutex.unlock();
    return fam;
}

// Load REG file from stream
reg_t load_reg(std::istream &in)
{
    reg_t reg;
    io_mutex.lock();
    for (in >> std::ws; !in.eof(); in >> std::ws)
        reg.push_back(load_refsam(in));
    io_mutex.unlock();
    return reg;
}

void set_encoding(unsigned char *data, size_t record_size, const size_t major, const size_t minor, const unsigned char g)
{
    assert(g <= 3); // TODO compiler directives for asserts
    const size_t index = major * record_size + minor / 4;
    const unsigned char value = data[index];
    const unsigned char mask = ~(3 << 2 * (minor % 4));
    data[index] = (g << 2 * (minor % 4)) | (value & mask);
}

void plink_t::set_haplotype(size_t indiv_index, size_t snp_index, int parent, bool bit)
{
    unsigned char old_g = this->genotype(indiv_index, snp_index);
    unsigned char new_g;
    if (bit)
        new_g = old_g | (1 << parent);
    else
        new_g = old_g & (~(1 << parent));
    if (new_g != old_g)
        this->set_genotype(indiv_index, snp_index, new_g);
}

bed_t *allocate_bed(size_t n_indiv, size_t n_snps, int orientation, bool phased)
{
    const size_t major = orientation == SNP_MAJOR_ORDER ? n_snps : n_indiv;
    const size_t minor = orientation == SNP_MAJOR_ORDER ? n_indiv : n_snps;
    const size_t record_size = bed_record_size(minor);
    const size_t data_size = major * record_size;
    return new bed_t(
               phased ? HAP_MAGIC1 : BED_MAGIC1,
               phased ? HAP_MAGIC2 : BED_MAGIC2,
               orientation,
               record_size,
               major * record_size,
               new unsigned char[data_size]);
}

bed_t *load_bed(FILE *fin, size_t n, size_t m, unsigned char *data)
{
    const int h1 = fgetc(fin);
    const int h2 = fgetc(fin);
    const int smo = fgetc(fin);
    if (feof(fin))
    {
        plink_exception_t ex;
        ex << "EOF while reading header";
        throw ex;
    }
    const size_t record_size = smo ? bed_record_size(n) : bed_record_size(m);
    const size_t expected_size = smo ? (m * record_size) : (n * record_size);
    if (data == NULL)
    {
        data = new unsigned char[expected_size];
        if (!data)
        {
            plink_exception_t ex;
            ex << "Failed to allocate data for binary ped";
            throw ex;
        }
    }
    const size_t bytesread = fread(data, sizeof(unsigned char), expected_size, fin);
    if (bytesread != expected_size)
    {
        plink_exception_t ex;
        ex << "Unexpected EOF reading binary data";
        delete[] data;
        throw ex;
    }
    fgetc(fin);
    if (!feof(fin))
    {
        plink_exception_t ex;
        ex << "Binary ped file larger than expected";
        delete[] data;
        throw ex;
    }
    return new bed_t(h1, h2, smo, record_size, bytesread, data);
}

size_t dump_bed(FILE *fout, const bed_t *bed)
{
    fputc(bed->header1, fout);
    fputc(bed->header2, fout);
    fputc(bed->orientation, fout);
    return fwrite(bed->data, sizeof(unsigned char), bed->size, fout);
}

// Index BIM file
std::map<std::string, size_t> bimindex(const bim_t &bim)
{
    std::map<std::string, size_t> bx;
    for (size_t m = 0; m < bim.size(); m++)
        bx[bim[m].rs] = m;
    return bx;
}

// Index FAM file
std::map<std::string, size_t> famindex(const fam_t &fam)
{
    std::map<std::string, size_t> fx;
    for (size_t i = 0; i < fam.size(); ++i)
        fx[fam[i].iid] = i;
    return fx;
}

std::map<std::string, std::string> regindex(const reg_t &reg)
{
    std::map<std::string, std::string> sm;
    for (size_t i = 0; i < reg.size(); ++i)
        sm[reg[i].iid] = reg[i].reg;
    return sm;
}

// Load(ped) to delete heap-allocated data
void clear_genotypes(std::vector<unsigned char *> &genotypes)
{
    for (size_t i = 0; i < genotypes.size(); ++i)
        delete[] genotypes[i];
}

std::ostream &operator<<(std::ostream &out, const plink_t &dataset)
{
    out << dataset.fam.size() << " individual(s), " << dataset.bim.size() << " SNP(s), ";
    if (dataset.bed)
    {
        out << "phased interpretation, " << (dataset.bed->orientation ? "SNP" : "individual") << "-major order";
    }
    else
        out << "(no bed loaded)";
    return out;
}

void dump(const plink_t &dataset, const std::string &bedfile, const std::string &bimfile, const std::string &famfile)
{
    plink_exception_t ex;
    std::ofstream famfout(famfile.c_str());
    if (!famfout)
    {
        ex << "failed to open .fam file " << famfile << " for writing";
        throw ex;
    }
    dump_fam(famfout, dataset.fam);
    famfout.close();
    std::ofstream bimfout(bimfile.c_str());
    if (!bimfout)
    {
        ex << "failed to open .bim file " << bimfile << " for writing";
        throw ex;
    }
    dump_bim(bimfout, dataset.bim);
    bimfout.close();
    FILE *bedfout = fopen(bedfile.c_str(), "w");
    if (!bedfout)
    {
        ex << "failed to open .bed file " << bedfile << " for writing";
        throw ex;
    }
    size_t bytes_written = dump_bed(bedfout, dataset.bed);
    if (bytes_written != dataset.bed->size)
    {
        ex << "failed to write .bed file " << bedfile;
        throw ex;
    }
    fclose(bedfout);
    return;
}

plink_t load(const std::string &stub)
{
    plink_t dataset;
    std::string famfile = stub + ".fam";
    std::string bimfile = stub + ".bim";
    std::string hapfile = stub + ".hap";

    // Load FAM file
    std::ifstream famfin(famfile.c_str());
    if (!famfin)
    {
        plink_exception_t ex;
        ex << "Failed to open .fam input file " << famfile;
        throw ex;
    }
    dataset.fam = load_fam(famfin);
    famfin.close();

    // Load BIM file
    std::ifstream bimfin(bimfile.c_str());
    if (!bimfin)
    {
        plink_exception_t ex;
        ex << "Failed to open .bim input file " << bimfile;
        throw ex;
    }
    dataset.bim = load_bim(bimfin);
    bimfin.close();

    // Load HAP file
    FILE *hapfin = fopen(hapfile.c_str(), "r");
    if (!hapfin)
    {
        plink_exception_t ex;
        ex << "Failed to open .hap input file " << hapfile;
        throw ex;
    }
    try
    {
        dataset.bed = load_bed(hapfin, dataset.fam.size(), dataset.bim.size());
        fclose(hapfin);
    }
    catch (plink_exception_t &ex)
    {
        plink_exception_t ex2;
        ex2 << ex.what() << " (in " << hapfile << ")";
        throw ex2;
    }
    if ((dataset.bed->header1 != BED_MAGIC1 && dataset.bed->header1 != HAP_MAGIC1) || (dataset.bed->header2 != BED_MAGIC2 && dataset.bed->header2 != HAP_MAGIC2) || (dataset.bed->orientation != INDIVIDUAL_MAJOR_ORDER && dataset.bed->orientation != SNP_MAJOR_ORDER))
    {
        plink_exception_t ex;
        ex << "incorrect header bytes in " << hapfile;
        throw ex;
    }
    return dataset;
}

} //namespace
