#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <plink.h>

/* Function prototypes */
int print_usage (const char *);

/* Externally defined variables */
extern int opterr, optopt, optind;

/* Main program function */

int
main (int argc, char *argv[])
{
    int c;
    int has_reg = 0;
    int is_phased = 0;
    char *instub;

    /* Get command line options */
    opterr = 0;
    while ((c = getopt (argc, argv, "rph")) != -1)
    {
        switch (c)
        {
            case 'r':
                has_reg = 1;
                break;
            case 'p':
                is_phased = 1;
                break;
            case 'h':
                print_usage(NULL);
                return EXIT_SUCCESS;
            case '?':
                fprintf(stderr, "[ERROR]: unknown option \"-%c\".\n", optopt);
                return EXIT_SUCCESS;
            default:
                return EXIT_SUCCESS;
        }
    }

    if (optind != argc - 1)
    {
        print_usage("Need PLINK input stub as mandatory argument");
        return EXIT_SUCCESS;
    }
    else
        instub = strdup(argv[optind]);

    /* Read in plink data set */
    plink_t *p = read_plink(instub, has_reg, is_phased);
    if (!p)
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
            char *hs0 = hap2str(p, i, 0);
            char *hs1 = hap2str(p, i, 1);
            printf("%s\t%s\t%s\n", p->fam[i].iid, res, hs0);
            printf("%s\t%s\t%s\n", p->fam[i].iid, res, hs1);
            free(hs0);
            free(hs1);
        }
    }

    return EXIT_SUCCESS;
}

int
print_usage (const char *msg)
{
    puts("Usage: test-api [OPTION]... [PLINK STUB]");
    puts("Test of libplink API\n");
    if (msg)
        printf("ERROR: %s\n\n", msg);
    puts("Options:");
    puts("  -p     Input data are phased");
    puts("  -r     Input PLINK stub includes a REG file");
    puts("  -h     Display this help message and exit");
    putchar('\n');  

    return 0;
}