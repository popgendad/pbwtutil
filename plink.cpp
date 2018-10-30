#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cstdlib>
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
