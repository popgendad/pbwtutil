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
read_plink (const char *instub, const int has_reg, const int is_phased)
{
    char *bimfile;
    char *famfile;
    char *bedfile;
    char *regfile;
    size_t regdum;
    size_t nsam;
    size_t nsnp;
    size_t stublen;
    plink_t *p;

    /* Define length of input data set stub string */
    stublen = strlen (instub);

    /* Allocate heap memory for plink data set */
    p = (plink_t *) malloc (sizeof(plink_t));
    if (p == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Allocate heap memory for bim file name string */
    bimfile = (char *) malloc (stublen + 4);
    if (bimfile == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Allocate heap memory for fam file name string */
    famfile = (char *) malloc (stublen + 4);
    if (famfile == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Allocate heap memory for bed/hap file name string */
    bedfile = (char *) malloc (stublen + 4);
    if (bedfile == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Allocate heap memory for reg file name string */
    regfile = (char *) malloc (stublen + 4);
    if (regfile == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Construct input file names from input stub */
    strcpy (bimfile, instub);
    strcpy (famfile, instub);
    strcpy (bedfile, instub);
    strcpy (regfile, instub);
    strcat (bimfile, ".bim");
    strcat (famfile, ".fam");
    strcat (regfile, ".reg");
    if (is_phased == TRUE)
        strcat (bedfile, ".hap");
    else
        strcat (bedfile, ".bed");

    /* Check if all input plink files exist */
    struct stat sinfo;
 
    if (stat (bimfile, &sinfo) != 0)
    {
        fprintf (stderr, "libplink [ERROR]: cannot access bimfile %s\n", bimfile);
        return NULL;
    }
    if (stat (famfile, &sinfo) != 0)
    {
        fprintf (stderr, "libplink [ERROR]: cannot access famfile %s\n", famfile);
        return NULL;
    }
    if (has_reg && stat (regfile, &sinfo) != 0)
    {
        fprintf (stderr, "libplink [ERROR]: cannot access regfile %s\n", regfile);
        return NULL;
    }
    if (stat (bedfile, &sinfo) != 0)
    {
        fprintf (stderr, "libplink [ERROR]: cannot access bed/hap file %s\n", bedfile);
        return NULL;
    }

    /* Read marker data from bim file */
    p->bim = read_bim (bimfile, &nsnp);
    
    /* Index the bim data set */
    p->bim_index = index_bim (p->bim, nsnp);
 
    /* Read the sample data from the fam file */
    p->fam = read_fam (famfile, &nsam);

    /* Index the fam data set */
    p->fam_index = index_fam (p->fam, nsam);

    /* Read the population data from the reg file */
    if (has_reg == TRUE)
    {
        p->reg = read_reg (regfile, &regdum);
        p->reg_index = index_reg (p->reg, regdum);
    }

    /* Read the genotype data from the bed/hap file */
    p->bed = read_bed (bedfile, nsam, nsnp, NULL);

    /* Assign metadata variables */
    p->nsam = nsam;
    p->nsnp = nsnp;

    /* Free heap memory for filename strings */
    free (bimfile);
    free (famfile);
    free (bedfile);
    free (regfile);

    return p;
}


void destroy_plink (plink_t *p)
{
    if (p->bed != NULL)
        destroy_bed (p->bed);
    if (p->bim != NULL)
        destroy_bim (p->bim, p->nsnp);
    if (p->fam != NULL)
        destroy_fam (p->fam, p->nsam);
    if (p->reg != NULL)
        destroy_reg (p->reg, p->nsam);
    if (p->bim_index != NULL)
        kh_destroy(integer, p->bim_index);
    if (p->fam_index != NULL)
        kh_destroy(integer, p->fam_index);
    if (p->reg_index != NULL)
        kh_destroy(integer, p->reg_index);
    free (p);

    return;
}

unsigned char *
hap2uchar (plink_t *p, const uint64_t i, const int parent)
{
	size_t j;
	unsigned char *str;

    /* Allocate heap memory for data array */   
    str = (unsigned char *) malloc (p->nsnp * sizeof(unsigned char));
    if (str == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Iterate through SNP positions */
	for (j = 0; j < p->nsnp; ++j)
        str[j] = (plink_haplotype (p->bed, i, j, parent)) ? '1' : '0';

	return str;
}


char *
hap2str (plink_t *p, const uint64_t i, const int parent)
{
    size_t j;
    char *str;

    /* Allocate heap memory for data string */
    str = (char *) malloc (p->nsnp + 1);
    if (str == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Iterate through SNP positions */
    for (j = 0; j < p->nsnp; ++j)
        str[j] = (plink_haplotype (p->bed, i, j, parent)) ? '1' : '0';

    /* Add null-terminating character at end of string */
    str[j] = '\0';

    return str;
}


char *
query_reg (plink_t *p, const char *iid)
{
    uint64_t i;
    char *result;
    khint_t k;

    /* Query reg hash for sample identifier */
    k = kh_get(integer, p->reg_index, iid);

    /* Store resulting reg entry */
    if (k != kh_end(p->reg_index))
    {
        i = kh_value(p->reg_index, k);
        result = strdup (p->reg[i].reg);
        return result;
    }
    else
        return NULL;
}


char *
query_pop (plink_t *p, const char *iid)
{
    uint64_t i;
    char *result;
    khint_t k;

    /* Query reg hash for sample identifier */
    k = kh_get(integer, p->reg_index, iid);

    /* Store resulting reg entry */
    if (k != kh_end(p->reg_index))
    {
        i = kh_value(p->reg_index, k);
        result = strdup (p->reg[i].pop);
        return result;
    }
    else
        return NULL;
}

uint64_t *
hap2ulong (plink_t *p, const uint64_t i, const int parent)
{
    size_t j;
    size_t nints;
    uint64_t *binarray;

    /* Define local variables */
    nints = p->nsnp / 64 + 1;
    binarray = (uint64_t *) malloc (nints * sizeof(uint64_t));

    /* Set all memory in binarray to zero */
    memset (binarray, 0, nints * sizeof(uint64_t));

    /* Set bits in binarray */
    for (j = 0; j < p->nsnp; ++j)
    {
        size_t k = j / 64;
        size_t m = j % 64;
        if (plink_haplotype (p->bed, i, j, parent) == 1)
            binarray[k] |= 1 << m;
    }

    return binarray;
}


bed_t *
read_bed (const char *bedfile, uint64_t nsam, uint64_t nsnp, unsigned char *data)
{
    bed_t *bed;
    FILE *fin;

    /* Allocate heap memory for bed_t structure */
    bed = (bed_t *) malloc (sizeof(bed_t));
    if (bed == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Open binary input stream */
    fin = fopen (bedfile, "rb");
    if (fin == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Read first three bytes */
    const int h1 = fgetc (fin);
    const int h2 = fgetc (fin);
    const int smo = fgetc (fin);

    /* No data */
    if (feof (fin) != 0)
    {
        fputs ("libplink [ERROR]: input bed/hap appears to be truncated\n", stderr);
        return NULL;
    }

    /* Determine memory block size based on orientation */
    const size_t record_size = smo ? bed_record_size(nsam) : bed_record_size(nsnp);
    const size_t expected_size = smo ? (nsnp * record_size) : (nsam * record_size);

    /* Allocate memory for data blob */
    if (data == NULL)
    {
        data = (unsigned char *) malloc (expected_size * sizeof(unsigned char));
        if (data == NULL)
            return NULL;
    }

    /* Read data blob into memory */
    const size_t bytesread = fread (data, sizeof(unsigned char), expected_size, fin);

    /* Possibly truncated file */
    if (bytesread != expected_size)
    {
        fputs ("Unexpected EOF reading binary data", stderr);
        free (data);
        return NULL;
    }

    /* Check if we are at the end of the file */
    fgetc (fin);
    if (feof (fin) == 0)
    {
        fputs ("Binary ped file larger than expected", stderr);
        free (data);
        return NULL;
    }

    /* Populate bed data structure */
    bed->header1 = h1;
    bed->header2 = h2;
    bed->orientation = smo;
    bed->record_size = record_size;
    bed->size = bytesread;
    bed->data = data;
    bed->phased = (bed->header2 == HAP_MAGIC2) ? TRUE : FALSE;

    return bed;
}


void destroy_bed (bed_t *bed)
{
    if (bed != NULL)
    {
        if (bed->data != NULL)
            free (bed->data);
        free (bed);
    }

    return;
}


bim_t *
read_bim (const char *bimfile, size_t *nl)
{
    size_t lines_read;
    size_t total_alloc = CHUNK_SIZE;
    char line[LINE_LENGTH];
    const char delim[] = " \t\n";
    bim_t *bim;
    FILE *fin;

    /* Open input file stream */
    fin = fopen (bimfile, "r");
    if (fin == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Initialize line counter */
    lines_read = 0;

    /* Allocate heap memory for bim_t structure */
    bim = (bim_t *) malloc (total_alloc * sizeof(bim_t));
    if (bim == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Iterate through lines of bim file */
    while (fgets (line, LINE_LENGTH, fin) != NULL)
    {
        ++lines_read;
        if (lines_read > total_alloc)
        {
            total_alloc += CHUNK_SIZE;
            bim = (bim_t *) realloc (bim, total_alloc * sizeof(bim_t));
        }

        /* Parse the bim entry */
        char *token;
        size_t i;

        /* Define index position */
        i = lines_read - 1;
   
        /* Parse tokens */
        token = strtok (line, delim);
        if (strcmp (token, "X") == 0)
            bim[i].ch = 23;
        else if (strcmp (token, "Y") == 0)
            bim[i].ch = 24;
        else if (strcmp (token, "XY") == 0)
            bim[i].ch = 25;
        else if (strcmp (token, "MT") == 0)
            bim[i].ch = 26;
        else
            bim[i].ch = atoi (token);
        token = strtok (NULL, delim);
        bim[i].rsid = strdup (token);
        token = strtok (NULL, delim);
        bim[i].cM = strtod (token, NULL);
        token = strtok (NULL, delim);
        bim[i].bp = strtoul (token, NULL, 0);
        token = strtok (NULL, delim);
        bim[i].a0 = token[0];
        token = strtok (NULL, delim);
        bim[i].a1 = token[0];
    }

    /* Resize bim array heap memory */
    bim = (bim_t *) realloc (bim, lines_read * sizeof(bim_t));

    /* Close the input file stream */
    fclose (fin);

    /* Set total number of lines read from bim file */
    *nl = lines_read;

    return bim;
}


void destroy_bim (bim_t *bim, const size_t nsnp)
{
    size_t i;

    if (bim != NULL)
    {
        for (i = 0; i < nsnp; ++i)
            if (bim[i].rsid != NULL)
                free (bim[i].rsid);
        free (bim);
    }

    return;
}


khash_t(integer) *
index_bim (const bim_t *bim, const size_t nl)
{
    int a;
    size_t i;
    khint_t k;
    khash_t(integer) *bx;

    /* Initialize bim hash */
    bx = kh_init(integer);

    /* Iterate through lines in bim and store indices keyed on rsid */
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
    size_t i;
    FILE *fout;

    /* Open output file stream */
    fout = fopen (outfile, "w");
    if (fout == NULL)
    {
        perror ("libplink [ERROR]");
        return -1;
    }

    for (i = 0; i < nl; ++i)
        fprintf (fout, "%d\t%s\t%lf\t%lu\t%c\t%c\n",
                 bim[i].ch, bim[i].rsid, bim[i].cM, bim[i].bp,
                 bim[i].a0, bim[i].a1);

    /* Close output file stream */
    fclose (fout);

    return 0;
}


fam_t *
read_fam (const char *famfile, size_t *nl)
{
    size_t lines_read;
    size_t total_alloc = CHUNK_SIZE;
    char line[LINE_LENGTH];
    const char delim[] = " \t\n";
    fam_t *fam;
    FILE *fin;

    /* Open fam input file stream */
    fin = fopen (famfile, "r");
    if (fin == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Allocate heap memory for fam_t structure */
    fam = (fam_t *) malloc (total_alloc * sizeof(fam_t));
    if (fam == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Initialize line counter */
    lines_read = 0;

    /* Iterate through the lines of the fam file */
    while (fgets (line, LINE_LENGTH, fin) != NULL)
    {
        ++lines_read;
        if (lines_read > total_alloc)
        {
            total_alloc += CHUNK_SIZE;
            fam = (fam_t *) realloc (fam, total_alloc * sizeof(fam_t));
        }

        /* Parse the fam entry */
        char *token;
        size_t i;

        /* Define index of position */
        i = lines_read - 1;

        /* Parse the tokens */
        token  = strtok (line, delim);
        fam[i].fid = strdup (token);
        token = strtok (NULL, delim);
        fam[i].iid = strdup (token);
        token = strtok (NULL, delim);
        fam[i].pid = strdup (token);
        token = strtok (NULL, delim);
        fam[i].mid = strdup (token);
        token = strtok (NULL, delim);
        fam[i].sex = strdup (token);
        token = strtok (NULL, delim);
        fam[i].phe = strdup (token);
    }

    /* Re-size the array holding fam data */
    fam = (fam_t *) realloc (fam, lines_read * sizeof(fam_t));

    /* Close the input file stream */
    fclose (fin);

    /* Assign the number of lines read from the fam file */
    *nl = lines_read;

    return fam;
}


void destroy_fam (fam_t *fam, const size_t nsam)
{
    size_t i;

    if (fam != NULL)
    {
        for (i = 0; i < nsam; ++i)
        {
            if (fam[i].fid != NULL)
                free (fam[i].fid);
            if (fam[i].iid != NULL)
                free (fam[i].iid);
            if (fam[i].pid != NULL)
                free (fam[i].pid);
            if (fam[i].mid != NULL)
                free (fam[i].mid);
            if (fam[i].sex != NULL)
                free (fam[i].sex);
            if (fam[i].phe != NULL)
                free (fam[i].phe);
        }
        free (fam);
    }

    return;
}


khash_t(integer) *
index_fam (const fam_t *fam, const size_t nl)
{
    int a;
    size_t i;
    khint_t k;
    khash_t(integer) *fx;

    /* Initialize fam hash */
    fx = kh_init(integer);

    /* Iterate through lines in bim and store indices keyed on sample identifier */
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
    size_t i;
    FILE *fout;

    fout = fopen (outfile, "w");
    if (fout == NULL)
    {
        perror ("libplink [ERROR]");
        return -1;
    }

    for (i = 0; i < nl; ++i)
        fprintf (fout, "%s\t%s\t%s\t%s\t%s\t%s\n",
                 fam[i].fid, fam[i].iid, fam[i].pid, fam[i].mid,
                 fam[i].sex, fam[i].phe);

    fclose (fout);

    return 0;
}


reg_t *
read_reg (const char *regfile, size_t *nl)
{
    size_t lines_read;
    size_t total_alloc = CHUNK_SIZE;
    reg_t *reg;
    FILE *fin;
    char line[LINE_LENGTH];
    const char delim[] = " \t\n";

    /* Open the input reg file stream */
    fin = fopen (regfile, "r");
    if (fin == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    reg = (reg_t *) malloc (total_alloc * sizeof(reg_t));
    if (reg == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Initialize the line counter */
    lines_read = 0;

    /* Iterate through the lines of the reg file */
    while (fgets (line, LINE_LENGTH, fin) != NULL)
    {
        ++lines_read;
        if (lines_read > total_alloc)
        {
            total_alloc += CHUNK_SIZE;
            reg = (reg_t *) realloc (reg, total_alloc * sizeof(reg_t));
        }

        /* Parse the reg entry */
        char *token;
        size_t i;

        /* Define the index of position */
        i = lines_read - 1;

        /* Parse the tokens */
        token = strtok (line, delim);
        reg[i].fid = strdup (token);
        token = strtok (NULL, delim);
        reg[i].iid = strdup (token);
        token = strtok (NULL, delim);
        reg[i].pid = strdup (token);
        token = strtok (NULL, delim);
        reg[i].mid = strdup (token);
        token = strtok (NULL, delim);
        reg[i].sex = strdup (token);
        token = strtok (NULL, delim);
        reg[i].phe = strdup (token);
        token = strtok (NULL, delim);
        reg[i].pop = strdup (token);
        token = strtok (NULL, delim);
        reg[i].reg = strdup (token);
    }

    /* Re-size the heap memory for the reg structure */
    reg = (reg_t *) realloc (reg, lines_read * sizeof(reg_t));

    /* Close the input file stream */
    fclose (fin);

    /* Assign the number of lines read from the reg file */
    *nl = lines_read;

    return reg;
}


void destroy_reg (reg_t *reg, const size_t nsam)
{
    size_t i;

    if (reg != NULL)
    {
        for (i = 0; i < nsam; ++i)
        {
            if (reg[i].fid != NULL)
                free (reg[i].fid);
            if (reg[i].iid != NULL)
                free (reg[i].iid);
            if (reg[i].pid != NULL)
                free (reg[i].pid);
            if (reg[i].mid != NULL)
                free (reg[i].mid);
            if (reg[i].sex != NULL)
                free (reg[i].sex);
            if (reg[i].phe != NULL)
                free (reg[i].phe);
            if (reg[i].pop != NULL)
                free (reg[i].pop);
            if (reg[i].reg != NULL)
                free (reg[i].reg);
        }
        free (reg);
    }

    return;
}


khash_t(integer) *
index_reg (const reg_t *reg, const size_t nl)
{
    int a;
    size_t i;
    khint_t k;
    khash_t(integer) *rx;

    rx = kh_init(integer);

    for (i = 0; i < nl; ++i)
    {
        k = kh_put(integer, rx, reg[i].iid, &a);
        kh_value(rx, k) = i;
    }

    return rx;
}


int
write_reg (const char *outfile, const reg_t *reg, const size_t nl)
{
    size_t i;
    FILE *fout;

    /* Open output file stream */
    fout = fopen (outfile, "w");
    if (fout == NULL)
    {
        perror ("libplink [ERROR]");
        return -1;
    }

    /* Print entry to file */
    for (i = 0; i < nl; ++i)
        fprintf (fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                 reg[i].fid, reg[i].iid, reg[i].pid, reg[i].mid,
                 reg[i].sex, reg[i].phe, reg[i].pop, reg[i].reg);

    /* Close the output file stream */
    fclose (fout);

    return 0;
}
