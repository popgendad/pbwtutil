#include "plink.h"


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
