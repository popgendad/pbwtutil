#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pbwtmaster.h"


int pbwt_convert_plink(cmd_t *c)
{
    int v = 0;
    size_t length = 0;
    char *outfile = NULL;
    pbwt_t *b = NULL;

    /* Construct outfile name */
    length = strlen(c->instub);
    outfile = (char *)malloc((length + 4) * sizeof(char));
    if (outfile == NULL)
    {
        fputs("pbwtmaster [ERROR]: memory allocation error\n", stderr);
        return -1;
    }
    strcpy(outfile, c->instub);
    strcat(outfile, ".pbwt");

    /* Import PBWT structure */
    b = pbwt_import_plink(c->instub);
    if (b == NULL)
    {
        fputs("pbwtmaster [ERROR]: problem importing PLINK stub\n", stderr);
        return -1;
    }

    /* Build the prefix and divergence arrays */
    v = pbwt_build(b);

    /* Write the pbwt to file */
    v = pbwt_write(outfile, b);
    if (v != 0)
    {
        fputs("pbwtmaster [ERROR]: Failed to write PBWT to disk\n", stderr);
        return -1;
    }

    /* Free memory for the data structure */
    pbwt_destroy(b);
    free(outfile);
    free(c->query);
    free(c);

    return 0;
}

int pbwt_convert_vcf(cmd_t *c)
{
    int v = 0;
    pbwt_t *b = NULL;

    /* Import PBWT structure */
    b = pbwt_import_vcf(c->instub, c->popmap);
    if (b == NULL)
    {
        fputs("pbwtmaster [ERROR]: problem importing VCF data\n", stderr);
        return -1;
    }

    /* Build the prefix and divergence arrays */
    v = pbwt_build(b);

    /* Write the pbwt to file */
    v = pbwt_write(c->outfile, b);
    if (v != 0)
    {
     	fputs("pbwtmaster [ERROR]: Failed to write PBWT to disk\n", stderr);
        return -1;
    }

    /* Clean up allocated memory */
    pbwt_destroy(b);
    free(c->query);
    free(c);

    return 0;
}
