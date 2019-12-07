#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/tbx.h>
#include <plink_lite.h>
#include "ancmatch.h"

khash_t(string) *read_popmap(const char *);
int check_popmap(const bcf_hdr_t *, const khash_t(string) *);

int pbwt_convert_plink(cmd_t *c)
{
    int v = 0;
    size_t i = 0;
    size_t len = 0;
    khint_t z = 0;
    size_t q = 0;
    char *outfile = NULL;
    plink_t *p = NULL;
    pbwt_t *b = NULL;

    /* Construct outfile name */
    i = strlen(c->instub);
    outfile = (char *)malloc(i + 4);
    strcpy(outfile, c->instub);
    strcat(outfile, ".pbwt");

    /* Initialize plink data structure */
    p = plink_init(c->instub, c->has_reg, c->is_phased);
    if (p == NULL)
    {
        return 1;
    }

    /* Initialize pbwt structure */
    b = pbwt_init(p->nsnp, 2 * p->nsam);
    if (b == NULL)
    {
        return 1;
    }

    /* Iterate through all samples in the fam/reg */
    for (i = 0; i < p->nsam; ++i)
    {
        memcpy(&b->data[TWODCORD(2*i, b->nsite, 0)],
            hap2uchar(p, i, 0),
            b->nsite * sizeof(unsigned char));
        memcpy(&b->data[TWODCORD(2*i+1, b->nsite, 0)],
            hap2uchar(p, i, 1),
            b->nsite * sizeof(unsigned char));
        len = strlen(p->fam[i].iid);
        b->sid[2*i] = (char *)malloc(len + 3);
        strcpy(b->sid[2*i], p->fam[i].iid);
        strcat(b->sid[2*i], ".1");
        b->sid[2*i+1] = (char *)malloc(len + 3);
        strcpy(b->sid[2*i+1], p->fam[i].iid);
        strcat(b->sid[2*i+1], ".2");
        z = kh_get(integer, p->reg_index, p->fam[i].iid);
        q = kh_value(p->reg_index, z);
        b->reg[2*i] = strdup(p->reg[q].reg);
        b->reg[2*i+1] = strdup(p->reg[q].reg);
    }

    /* Iterate through all sites in the bim */
    for (i = 0; i < p->nsnp; ++i)
    {
        b->rsid[i] = strdup (p->bim[i].rsid);
        b->cm[i] = p->bim[i].cM;
        b->chr[i] = p->bim[i].ch;
    }

    /* Build the prefix and divergence arrays */
    v = pbwt_build(b);

    /* Write the pbwt to file */
    v = pbwt_write(outfile, b);
    if (v != 0)
    {
        fputs("ancmatch [ERROR]: Failed to write pbwt to disk", stderr);
        return 1;
    }

    /* Free memory for the data structure */
    pbwt_destroy(b);

    return 0;
}

int pbwt_convert_vcf(cmd_t *c)
{
    int32_t i = 0;
    int32_t j = 0;
    int v = 0;
    int32_t nseq = 0;
    uint64_t nsites = 0;
    int32_t nsam;
    const char **seq;
    tbx_t *tbx = NULL;
    bcf_hdr_t *hdr = NULL;
    htsFile *infile = NULL;
    hts_idx_t *idx = NULL;
    khash_t(string) *popdb = NULL;
    bcf1_t *rec = NULL;
    pbwt_t *b = NULL;

    /* Read the popmap file into a hash table */
    if ((popdb = read_popmap(c->popmap)) == NULL)
    {
        return 1;
    }

    /* Open the input file stream */
    if ((infile = hts_open(c->instub, "r")) == NULL)
    {
        return 1;
    }

    /* Read the VCF header into memory */
    if ((hdr = bcf_hdr_read(infile)) == NULL)
    {
        return 1;
    }

    /* Load VCF index file */
    if (hts_get_format(infile)->format == vcf)
    {
        if ((tbx = tbx_index_load(c->instub)) == NULL)
        {
            return 1;
        }
    }

    /* Get number of sites in VCF */
    seq = tbx ? tbx_seqnames(tbx, &nseq) : bcf_index_seqnames(idx, hdr, &nseq);
    for (i = 0; i < nseq; ++i)
    {
        uint64_t records = 0;
        uint64_t v = 0;
        hts_idx_get_stat(tbx ? tbx->idx : idx, i, &records, &v);
        nsites += records;
    }

    /* Check whether all VCF input samples are in popmap file */
    if ((nsam = check_popmap(hdr, popdb)) < 0)
    {
        return 1;
    }

    /* Initialize memory to hold single VCF record */
    if ((rec = bcf_init()) == NULL)
    {
        bcf_hdr_destroy(hdr);
        hts_close(infile);
        return 1;
    }

    /* Initialize pbwt structure */
    b = pbwt_init(nsites, 2 * nsam);
    if (b == NULL)
    {
        return 1;
    }

    size_t site_counter = 0;

    /* Read through single site entry in VCF */
    while (bcf_read(infile, hdr, rec) == 0)
    {
        const char *chr = NULL;
        int32_t ngt = 0;
        int32_t *gt = NULL;
        khint_t it = 0;
        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_INFO);
        chr = bcf_seqname(hdr, rec);
        b->chr[site_counter] = atoi(chr);
        b->rsid[site_counter] = strdup(rec->d.id);
        b->cm[site_counter] = rec->d.info[0].v1.f;
        ngt = bcf_get_genotypes(hdr, rec, &gt, &ngt) / nsam;

        /* Iterate through sample genotypes */
        for (i = 0; i < nsam; ++i)
        {
            const char *sid = hdr->samples[i];
            it = kh_get(string, popdb, sid);
            const char *pop = kh_value(popdb, it);
            int32_t *ptr = gt + i * ngt;
            const size_t len = strlen(sid);
            b->sid[2*i] = (char *)malloc(len + 3);
            b->sid[2*i+1] = (char *)malloc(len + 3);
            strcpy(b->sid[2*i], sid);
            strcat(b->sid[2*i], ".1");
            strcpy(b->sid[2*i+1], sid);
            strcat(b->sid[2*i+1], ".2");
            b->reg[2*i] = strdup(pop);
            b->reg[2*i+1] = strdup(pop);
            for (j = 0; j < ngt; ++j)
            {
                if (bcf_gt_is_missing(ptr[j]))
                {
                    break;
                }

                if (ptr[j] == bcf_int32_vector_end)
                {
                    break;
                }

                if (j % 2 == 0)
                {
                    if (bcf_gt_allele(ptr[j]) == 1)
                        b->data[TWODCORD(2*i, b->nsite, site_counter)] = '1';
                    else
                        b->data[TWODCORD(2*i, b->nsite, site_counter)] = '0';
                }
                else
                {
                    if (bcf_gt_allele(ptr[j]) == 1)
                        b->data[TWODCORD(2*i+1, b->nsite, site_counter)] = '1';
                    else
                        b->data[TWODCORD(2*i+1, b->nsite, site_counter)] = '0';
                }
            }
        }
        site_counter += 1;
    }

    /* Build the prefix and divergence arrays */
    v = pbwt_build(b);

    /* Write the pbwt to file */
    v = pbwt_write(c->outfile, b);
    if (v != 0)
    {
     	fputs("ancmatch [ERROR]: Failed to write pbwt to disk", stderr);
        return 1;
    }

    /* Free memory for the data structure */
    pbwt_destroy(b);

    /* Close up shop */
    free(seq);
    bcf_hdr_destroy(hdr);
    hts_close(infile);
    bcf_destroy(rec);
    tbx_destroy(tbx);
    hts_idx_destroy(idx);

    return 0;
}


int check_popmap(const bcf_hdr_t *h, const khash_t(string) *pdb)
{
    int k = 0;
    khint_t it = 0;
    const int32_t nsam = bcf_hdr_nsamples(h);

    /* Check the sample identifiers in the VCF header against the database */
    for (k = 0; k < nsam; ++k)
    {
        const char *sid = h->samples[k];
        it = kh_get(string, pdb, sid);
        if (it == kh_end(pdb))
        {
            fprintf(stderr, "ancmatch [ERROR]: sample not in population database: %s\n", sid);
            return -1;
        }
    }
    return nsam;
}


khash_t(string) *read_popmap(const char *popfile)
{
    khash_t(string) *popdb = NULL;
    khint_t it = 0;
    char sid[256];
    char pop[256];
    int counter = 0;
    FILE *instream;

    /* Initialize sample->population hash table */
    popdb = kh_init(string);

    /* Open popmap input file stream */
    instream = fopen(popfile, "r");
    if (!instream)
    {
        return NULL;
    }

    while (!feof(instream))
    {
        int absent = 0;
        int ns = 0;

        ns = fscanf(instream, "%s\t%s\n", sid, pop);
        if (ns != 2)
        {
            fputs("ancmatch [ERROR] cannot read popmap file\n", stderr);
            return NULL;
        }
        it = kh_put(string, popdb, sid, &absent);
        if (absent)
        {
            kh_key(popdb, it) = strdup(sid);
        }
        kh_value(popdb, it) = strdup(pop);
    	counter++;
    }

    fclose(instream);
    fprintf(stderr, "ancmatch [INFO]: %d entries from %s read into pop database\n", counter, popfile);

    return popdb;
}
