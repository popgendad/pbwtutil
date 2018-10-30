#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "plink.h"

extern int opterr, optopt, optind;

/* ./test-api -p -r EV62-phase.3881 */

int
main (int argc, char *argv[])
{
    int c;
    int has_reg = 0;
    int is_phased = 0;
    char *instub;

    /* Get command line options */
    opterr = 0;
    while ((c = getopt (argc, argv, "rp")) != -1)
    {
        switch (c)
        {
            case 'r':
                has_reg = 1;
                break;
            case 'p':
                is_phased = 1;
                break;
            case '?':
                fprintf(stderr, "ERROR: unknown option \"-%c\".\n", optopt);
            default:
                return EXIT_SUCCESS;
        }
    }
    instub = strdup(argv[optind]);


    /* Read in plink data set */
    plink_t *p = read_plink(instub, has_reg, is_phased);
    if (p == NULL)
        return EXIT_FAILURE;


    /* Example: How to use query_reg() function
     * to read data from test plink set and print
     * out all haplotypes in Beringia reference panel */
    size_t i = 0;
    char qreg[] = "Beringia";

    /* Iterate through all samples in the fam/reg */
    for (i = 0; i < p->nsam; ++i)
    {
        char *res = query_reg(p->reg, p->fam[i].iid);
        if ((res != NULL) && (strcmp(qreg , res) == 0))
        {
            printf("%s\t%s\t", p->fam[i].iid, res);
            char *hs0 = hap2str(p, i, 0);
            printf("%s\n", hs0);
            printf("%s\t%s\t", p->fam[i].iid, res);
            char *hs1 = hap2str(p, i, 1);
            printf("%s\n", hs1);
            free(hs0);
            free(hs1);
        }
    }

    return EXIT_SUCCESS;
}
