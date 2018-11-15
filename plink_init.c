#include <sys/stat.h>
#include <sys/types.h>
#include "plink.h"

plink_t *
plink_init (const char *instub, const int has_reg, const int is_phased)
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
    p->bim = bim_read (bimfile, &nsnp);
    if (p->bim == NULL)
    {
        fputs ("libplink [ERROR]: problem reading bim file", stderr);
        return NULL;
    }
    
    /* Index the bim data set */
    p->bim_index = bim_index (p->bim, nsnp);
    if (p->bim_index == NULL)
    {
        fputs ("libplink [ERROR]: problem indexing bim file", stderr);
        return NULL;
    }
 
    /* Read the sample data from the fam file */
    p->fam = fam_read (famfile, &nsam);
    if (p->fam == NULL)
    {
        fputs ("libplink [ERROR]: problem reading fam file", stderr);
        return NULL;
    }

    /* Index the fam data set */
    p->fam_index = fam_index (p->fam, nsam);
    if (p->fam_index == NULL)
    {
        fputs ("libplink [ERROR]: problem indexing fam file", stderr);
        return NULL;
    }

    /* Read the population data from the reg file */
    if (has_reg == TRUE)
    {
        p->reg = reg_read (regfile, &regdum);
        if (p->reg == NULL)
        {
            fputs ("libplink [ERROR]: problem reading reg file", stderr);
            return NULL;
        }
        p->reg_index = reg_index (p->reg, regdum);
        if (p->reg_index == NULL)
        {
            fputs ("libplink [ERROR]: problem indexing reg file", stderr);
            return NULL;
        }
    }

    /* Read the genotype data from the bed/hap file */
    p->bed = bed_read (bedfile, nsam, nsnp, NULL);
    if (p->bed == NULL)
    {
        fputs ("libplink [ERROR]: problem reading bed/hap file", stderr);
        return NULL;
    }

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
