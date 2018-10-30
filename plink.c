#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#include "plink.h"

#define CHUNK_SIZE 1000
#define LINE_LENGTH 1024

/* Unknown allele value in .bim file */
const char UNKNOWN_VARIANT = '0';

plink_t *
read_plink(const char *instub, int phased)
{
    plink_t *p = malloc(sizeof(plink_t));
    size_t regdum, nsam, nsnp;
    size_t stublen = strlen(instub);
    char *bimfile = malloc(stublen + 4);
    char *famfile = malloc(stublen + 4);
    char *bedfile = malloc(stublen + 4);
    char *regfile = malloc(stublen + 4);

    /* Construct input file names from input stub */
    strcpy(bimfile, instub);
    strcpy(famfile, instub);
    strcpy(bedfile, instub);
    strcpy(regfile, instub);
 
    strcat(bimfile, ".bim");
    strcat(famfile, ".fam");
    strcat(regfile, ".reg");
    if (phased)
        strcat(bedfile, ".hap");
    else
        strcat(bedfile, ".bed");

    /* Check if all input plink files exist */
    struct stat sinfo;
    if (stat(bimfile, &sinfo) != 0)
    {
        fprintf(stderr, "ERROR: problem accessing bimfile %s\n", bimfile);
        return NULL;
    }
    if (stat(famfile, &sinfo) != 0)
    {
        fprintf(stderr, "ERROR: problem accessing famfile %s\n", famfile);
        return NULL;
    }
    if (stat(regfile, &sinfo) != 0)
    {
        fprintf(stderr, "ERROR: problem accessing regfile %s\n", regfile);
        return NULL;
    }
    if (stat(bedfile, &sinfo) != 0)
    {
        fprintf(stderr, "ERROR: problem accessing bedfile %s\n", bedfile);
        return NULL;
    }

    /* Read data from files */
    p->bim = read_bim(bimfile, &nsnp);
    p->fam = read_fam(famfile, &nsam);
    p->reg = read_reg(regfile, &regdum);
    p->bed = read_bed(bedfile, nsam, nsnp, NULL);
    p->nsam = nsam;
    p->nsnp = nsnp;

    /* Free heap memory for filename strings */
    free(bimfile);
    free(famfile);
    free(bedfile);
    free(regfile);

    return p;
}

bed_t *
read_bed (const char *bedfile, uint64_t nsam, uint64_t nsnp, unsigned char *data)
{
    bed_t *bed = malloc(sizeof(bed_t));
    FILE *fin;
    fin = fopen(bedfile, "r");
    const int h1 = fgetc(fin);
    const int h2 = fgetc(fin);
    const int smo = fgetc(fin);

    if (feof(fin))
        return 0;

    const size_t record_size = smo ? bed_record_size(nsam) : bed_record_size(nsnp);
    const size_t expected_size = smo ? (nsnp * record_size) : (nsam * record_size);

    if (data == NULL)
    {
        data = malloc(expected_size * sizeof(unsigned char));
        if (!data)
            return 0;
    }

    const size_t bytesread = fread(data, sizeof(unsigned char), expected_size, fin);
    if (bytesread != expected_size)
    {
        fprintf(stderr, "Unexpected EOF reading binary data");
        free(data);
        return 0;
    }

    fgetc(fin);
    if (!feof(fin))
    {
        fprintf(stderr, "Binary ped file larger than expected");
        free(data);
        return 0;
    }
 
    bed->header1 = h1;
    bed->header2 = h2;
    bed->orientation = smo;
    bed->record_size = record_size;
    bed->size = bytesread;
    bed->data = data;
    if (bed->header2 == HAP_MAGIC2)
        bed->phased = 1;
    else
        bed->phased = 0;
 
    return bed;
}

bim_t *
read_bim (const char *bimfile, size_t *nl)
{
    size_t lc = 0;
    size_t lalloc = CHUNK_SIZE;
    char line[LINE_LENGTH];
    const char delim[] = " \t";
    bim_t *bim;
    FILE *fin;

    fin = fopen(bimfile, "r");
    if (fin == NULL)
        return 0;
 
    bim = malloc(lalloc * sizeof(bim_t));
 
    while (fgets(line, LINE_LENGTH, fin) != NULL)
    {
        ++lc;
        if (lc > lalloc)
        {
            lalloc += CHUNK_SIZE;
            bim = realloc(bim, lalloc * sizeof(bim_t));
        }

        char *token = strtok(line, delim);
        if (strcmp(token, "X") == 0)
            bim[lc-1].ch = 23;
        else if (strcmp(token, "Y") == 0)
            bim[lc-1].ch = 24;
        else if (strcmp(token, "XY") == 0)
            bim[lc-1].ch = 25;
        else if (strcmp(token, "MT") == 0)
            bim[lc-1].ch = 26;
        else
            bim[lc-1].ch = atoi(token);
        token = strtok(NULL, delim);
        bim[lc-1].rsid = strdup(token);      
        token = strtok(NULL, delim);
        bim[lc-1].cM = strtod(token, NULL);
        token = strtok(NULL, delim);
        bim[lc-1].bp = strtoul(token, NULL, 0);
        token = strtok(NULL, delim);
        bim[lc-1].a0 = token[0];
        token = strtok(NULL, delim);
        bim[lc-1].a1 = token[0];
    }

    bim = realloc(bim, lalloc * sizeof(bim_t));

    fclose(fin);

#ifdef DEBUG
    printf("%lu SNPs read from bim file %s\n", lc, bimfile);
#endif

    /* Index the bim data set */
    bim->index = index_bim(bim, lc);

    *nl = lc;
    return bim;
}

khash_t(integer) *
index_bim (const bim_t *bim, const size_t nl)
{
    int a = 0;
    size_t i = 0;
    khint_t k = 0;
    khash_t(integer) *bx = NULL;

    bx = kh_init(integer);

    for (i = 0; i < nl; ++i)
    {
        k = kh_put(integer, bx, bim[i].rsid, &a);
        kh_value(bx, k) = i;
    }

    return bx;
}

int
write_bim (const char *outfile, const bim_t *bim, const size_t nl)
{
    size_t i = 0;
    FILE *fout;

    fout = fopen(outfile, "w");

    for (i = 0; i < nl; ++i)
        fprintf(fout, "%d\t%s\t%lf\t%lu\t%c\t%c\n",
                bim[i].ch, bim[i].rsid, bim[i].cM, bim[i].bp,
                bim[i].a0, bim[i].a1);

    fclose(fout);

    return 0;
}

fam_t *
read_fam (const char *famfile, size_t *nl)
{
    size_t lc = 0;
    size_t lalloc = CHUNK_SIZE;
    char line[LINE_LENGTH];
    const char delim[] = " \t";
    fam_t *fam;
    FILE *fin;

    fin = fopen(famfile, "r");
    if (fin == NULL)
        return 0;

    fam = malloc(lalloc * sizeof(fam_t));

    while (fgets(line, LINE_LENGTH, fin) != NULL)
    {
        ++lc;
        if (lc > lalloc)
        {
            lalloc += CHUNK_SIZE;
            fam = realloc(fam, lalloc * sizeof(fam_t));
        }

        char *token = strtok(line, delim);
        fam[lc-1].fid = strdup(token);
        token = strtok(NULL, delim);
        fam[lc-1].iid = strdup(token);     
        token = strtok(NULL, delim);
        fam[lc-1].pid = strdup(token);
        token = strtok(NULL, delim);
        fam[lc-1].mid = strdup(token);
        token = strtok(NULL, delim);
        fam[lc-1].sex = strdup(token);
        token = strtok(NULL, delim);
        fam[lc-1].phe = strdup(token);
    }

    fam = realloc(fam, lalloc * sizeof(fam_t));

    fclose(fin);

#ifdef DEBUG
    fprintf(stderr, "%lu diploid samples read from fam file %s\n", lc, famfile);
#endif

    /* Index the fam data set */
    fam->index = index_fam(fam, lc);

    *nl = lc;
    return fam;
}

khash_t(integer) *
index_fam (const fam_t *fam, const size_t nl)
{
    int a = 0;
    size_t i = 0;
    khint_t k = 0;
    khash_t(integer) *fx = NULL;

    fx = kh_init(integer);

    for (i = 0; i < nl; ++i)
    {
        k = kh_put(integer, fx, fam[i].iid, &a);
        kh_value(fx, k) = i;
    }

    return fx;   
}

int
write_fam (const char *outfile, const fam_t *fam, const size_t nl)
{
    size_t i = 0;
    FILE *fout;

    fout = fopen(outfile, "w");

    for (i = 0; i < nl; ++i)
        fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\n",
        fam[i].fid, fam[i].iid, fam[i].pid, fam[i].mid,
        fam[i].sex, fam[i].phe);

    fclose(fout);

    return 0;
}

reg_t *
read_reg (const char *regfile, size_t *nl)
{
    reg_t *reg;
    FILE *fin;
    char line[LINE_LENGTH];
    size_t lc = 0;
    size_t lalloc = CHUNK_SIZE;
    const char delim[] = " \t";

    fin = fopen(regfile, "r");
    if (fin == NULL)
        return 0;

    reg = malloc(lalloc * sizeof(reg_t));

    while (fgets(line, LINE_LENGTH, fin) != NULL)
    {
        ++lc;
        if (lc > lalloc)
        {
            lalloc += CHUNK_SIZE;
            reg = realloc(reg, lalloc * sizeof(reg_t));
        }

        char *token = strtok(line, delim);
        reg[lc-1].fid = strdup(token);
        token = strtok(NULL, delim);
        reg[lc-1].iid = strdup(token);     
        token = strtok(NULL, delim);
        reg[lc-1].pid = strdup(token);
        token = strtok(NULL, delim);
        reg[lc-1].mid = strdup(token);
        token = strtok(NULL, delim);
        reg[lc-1].sex = strdup(token);
        token = strtok(NULL, delim);
        reg[lc-1].phe = strdup(token);
        token = strtok(NULL, delim);
        reg[lc-1].pop = strdup(token);
        token = strtok(NULL, delim);
        reg[lc-1].reg = strdup(token);
    }

    reg = realloc(reg, lalloc * sizeof(reg_t));

    fclose(fin);

#ifdef DEBUG
    fprintf(stderr, "%lu samples read from reg file %s\n", lc, regfile);
#endif

    /* Index the reg data set */
    reg->index = index_reg(reg, lc);

    *nl = lc;
    return reg;
}

khash_t(integer) *
index_reg (const reg_t *reg, const size_t nl)
{
    int a = 0;
    size_t i = 0;
    khint_t k = 0;
    khash_t(integer) *rx = NULL;

    rx = kh_init(integer);

    for (i = 0; i < nl; ++i)
    {
        k = kh_put(integer, rx, reg[i].iid, &a);
        kh_value(rx, k) = i;
    }

    return rx;     
}

char *
query_reg(const char *iid, reg_t *r)
{
    uint64_t i;
    char *result;
    khint_t k = 0;
    k = kh_get(integer, r->index, iid);
    if (k != kh_end(r->index))
        i = kh_value(r->index, k);
    result = strdup(r[i].reg);

    return result;
}

int
write_reg (const char *outfile, const reg_t *reg, const size_t nl)
{
    size_t i = 0;
    FILE *fout;

    fout = fopen(outfile, "w");

    for (i = 0; i < nl; ++i)
        fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                reg[i].fid, reg[i].iid, reg[i].pid, reg[i].mid, 
                reg[i].sex, reg[i].phe, reg[i].pop, reg[i].reg);

    fclose(fout);

    return 0;
}
