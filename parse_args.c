#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "pbwtmaster.h"

/* Get version information */
#define XSTR(x) #x
#define STR(x) XSTR(x)
#ifdef VERSION
const char Version[] = STR(VERSION);
#endif

/* Externally defined variables */
extern int opterr, optopt, optind;
extern char *optarg;

/* Local function prototypes */
int parse_coancestry(int, char **, cmd_t *);
int parse_convert(int, char **, cmd_t *);
int parse_match(int, char **, cmd_t *);
int parse_pileup(int, char **, cmd_t *);
int parse_summary(int, char **, cmd_t *);
int parse_view(int, char **, cmd_t *);
int print_main_usage(const char *);
int print_coancestry_usage(const char *);
int print_convert_usage(const char *);
int print_match_usage(const char *);
int print_pileup_usage(const char *);
int print_summary_usage(const char *);
int print_view_usage(const char *);
void print_version(void);

cmd_t *parse_args(int argc, char *argv[])
{
    char *mode = NULL;
    cmd_t *c = NULL;
    int (*parse_func)(int, char **, cmd_t *);

    opterr = 0;

    /* Allocate memory for command line data structure */
    c = (cmd_t *)malloc(sizeof(cmd_t));
    if (c == NULL)
    {
        fputs("pbwtmaster [ERROR]: memory allocation failure\n", stderr);
        return NULL;
    }

    /* Initialize default values for command line options */
    c->has_reg = 0;
    c->is_phased = 0;
    c->with_vcf = 0;
    c->match_all = 0;
    c->nohaps = 0;
    c->minlen = 0.5;
    c->only_sites = 0;
    c->print_sites = 0;
    c->count_only = 0;
    c->reg_count = 0;
    c->adjlist = 0;
    c->set_match = 0;
    c->out_diploid = 0;

    /* Get mode argument */
    if (argv[1])
    {
        mode = strdup(argv[1]);
    }
    else
    {
        print_main_usage(NULL);
        return NULL;
    }

    /* Iterate command line variables */
    argc--;
    argv++;

    /* Determine run-time mode */
    if (strcmp(mode, "coancestry") == 0)
    {
        c->mode = COANCESTRY;
        c->mode_func = &pbwt_coancestry;
        parse_func = &parse_coancestry;
    }
    else if (strcmp(mode, "convert") == 0)
    {
        c->mode = CONVERT;
        parse_func = &parse_convert;
    }
    else if (strcmp(mode, "match") == 0)
    {
        c->mode = MATCH;
        c->mode_func = &pbwt_match;
        parse_func = &parse_match;
    }
    else if (strcmp(mode, "pileup") == 0)
    {
        c->mode = PILEUP;
        c->mode_func = &pbwt_pileup;
        parse_func = &parse_pileup;
    }
    else if (strcmp(mode, "summary") == 0)
    {
        c->mode = SUMMARY;
        c->mode_func = &pbwt_summary;
        parse_func = &parse_summary;
    }
    else if (strcmp(mode, "view") == 0)
    {
        c->mode = VIEW;
        c->mode_func = &pbwt_view;
        parse_func = &parse_view;
    }
    else
    {
        /* No match for mode string */
        if (mode)
        {
            char msg[100];
            sprintf(msg, "pbwtmaster [ERROR]: unknown command \'%s\'", mode);
            print_main_usage(msg);
        }
        else
        {
            print_main_usage("pbwtmaster [ERROR]: need to specify a command");
        }

        return NULL;
    }

    /* Call the command-specific parsing function */
    int r = (*parse_func)(argc, argv, c);
    if (r < 0)
    {
        return NULL;
    }

    if (strcmp(mode, "convert") == 0)
    {
        if (c->with_vcf == 0)
        {
            c->mode_func = &pbwt_convert_plink;
        }
        else
        {
            c->mode_func = &pbwt_convert_vcf;
        }
    }

    /* Remove mode string from heap */
    free(mode);

    return c;
}

int parse_coancestry(int argc, char *argv[], cmd_t *c)
{
    int g = 0;
    char msg[100];

    while (1)
    {
        int option_index = 0;

        /* Declare the option table */
        static struct option long_options[] =
        {
            { "diploid", no_argument,       NULL, 'd' },
            { "adjlist", no_argument,       NULL, 'a' },
            { "set",     no_argument,       NULL, 's' },
            { "sites",   no_argument,       NULL, 'p' },
            { "count",   no_argument,       NULL, 'c' },
            { "minlen",  required_argument, NULL, 'm' },
            { "version", no_argument,       NULL, 'v' },
            { "help",    no_argument,       NULL, 'h' },
            {0, 0, 0, 0}
        };

        /* Parse the option */
        g = getopt_long(argc, argv, "daspcvhm:", long_options, &option_index);

        /* We are at the end of the options */
        if (g == -1)
        {
            break;
        }

        /* Assign the option to variables */
        switch(g)
        {
            case 'd':
                c->out_diploid = 1;
                break;
            case 'a':
                c->adjlist = 1;
                break;
            case 'm':
                c->minlen = atof(optarg);
                break;
            case 's':
                c->set_match = 1;
                break;
            case 'c':
                c->count_only = 1;
                break;
            case 'p':
                c->print_sites = 1;
                break;
            case 'v':
                print_version();
                return -1;
            case 'h':
                print_coancestry_usage(NULL);
                return -1;
            case '?':
                sprintf(msg, "pbwtmaster [ERROR]: unknown option \"-%c\".\n", optopt);
                print_coancestry_usage(msg);
                return -1;
            default:
                print_coancestry_usage(NULL);
                return -1;
        }
    }

    /* Parse non-optioned arguments */
    if (optind != argc - 1)
    {
        print_coancestry_usage("pbwtmaster [ERROR]: need input file name as mandatory argument");
        return -1;
    }
    else
    {
        c->instub = strdup(argv[optind]);
    }

    return 0;
}

int parse_convert(int argc, char *argv[], cmd_t *c)
{
    int g = 0;
    char msg[100];

    while (1)
    {
        int option_index = 0;

        /* Declare the option table */
        static struct option long_options[] =
        {
            { "vcf",     no_argument,       NULL, 'c' },
            { "map",     required_argument, NULL, 'm' },
            { "reg",     no_argument,       NULL, 'r' },
            { "out",     required_argument, NULL, 'o' },
            { "phased",  no_argument,       NULL, 'p' },
            { "version", no_argument,       NULL, 'v' },
            { "help",    no_argument,       NULL, 'h' },
            {0, 0, 0, 0}
        };

        /* Parse the option */
        g = getopt_long(argc, argv, "vhcrpm:o:", long_options, &option_index);

        /* We are at the end of the options */
        if (g == -1)
        {
            break;
        }

        /* Assign the option to variables */
        switch(g)
        {
            case 'c':
                c->with_vcf = 1;
                break;
            case 'r':
                c->has_reg = 1;
                break;
            case 'o':
                c->outfile = strdup(optarg);
                break;
            case 'm':
                c->popmap = strdup(optarg);
                break;
            case 'p':
                c->is_phased = 1;
                break;
            case 'v':
                print_version();
                return -1;
            case 'h':
                print_convert_usage(NULL);
                return -1;
            case '?':
                sprintf(msg, "pbwtmaster [ERROR]: unknown option \"-%c\".\n", optopt);
                print_convert_usage(msg);
                return -1;
            default:
                print_convert_usage(NULL);
                return -1;
        }
    }

    /* Parse non-optioned arguments */
    if (optind != argc - 1)
    {
        print_convert_usage("pbwtmaster [ERROR]: need input stub as mandatory argument");
        return -1;
    }
    else
    {
        c->instub = strdup(argv[optind]);
    }

    /* --out switch is actually mandatory */
    if (!c->outfile)
    {
        print_convert_usage("pbwtmaster [ERROR]: --out <STR> is a mandatory argument for convert");
        return -1;
    }

    return 0;
}

int parse_match(int argc, char *argv[], cmd_t *c)
{
    int g = 0;
    char msg[100];

    while (1)
    {
        int option_index = 0;

        /* Declare the option table */
        static struct option long_options[] =
        {
            { "query",   required_argument, NULL, 'q' },
            { "minlen",  required_argument, NULL, 'm' },
            { "set",     no_argument,       NULL, 's' },
            { "sites",   no_argument,       NULL, 'p' },
            { "all",     no_argument,       NULL, 'a' },
            { "version", no_argument,       NULL, 'v' },
            { "help",    no_argument,       NULL, 'h' },
            {0, 0, 0, 0}
        };

        /* Parse the option */
        g = getopt_long(argc, argv, "vhapq:m:", long_options, &option_index);

        /* We are at the end of the options */
        if (g == -1)
        {
            break;
        }

        /* Assign the option to variables */
        switch(g)
        {
            case 'q':
                c->query = strdup(optarg);
                break;
            case 'p':
                c->print_sites = 1;
                break;
            case 'm':
                c->minlen = atof(optarg);
                break;
            case 'a':
                c->match_all = 1;
                break;
            case 's':
                c->set_match = 1;
                break;
            case 'v':
                print_version();
                return -1;
            case 'h':
                print_match_usage(NULL);
                return -1;
            case '?':
                sprintf(msg, "pbwtmaster [ERROR]: unknown option \"-%c\".\n", optopt);
                print_match_usage(msg);
                return -1;
            default:
                print_match_usage(NULL);
                return -1;
        }
    }

    /* Parse non-optioned arguments */
    if (optind != argc - 1)
    {
        print_match_usage("pbwtmaster [ERROR]: need input file name as mandatory argument");
        return -1;
    }
    else
    {
        c->instub = strdup(argv[optind]);
    }

    /* Check that a query sequence has been specified */
    if (c->query == NULL)
    {
        print_match_usage("pbwtmaster [ERROR]: --query option is mandatory");
        return -1;
    }

    return 0;
}

int parse_pileup(int argc, char *argv[], cmd_t *c)
{
    int g = 0;
    char msg[100];

    while (1)
    {
        int option_index = 0;

        /* Declare the option table */
        static struct option long_options[] =
        {
            { "query",   required_argument, NULL, 'q' },
            { "minlen",  required_argument, NULL, 'm' },
            { "set",     no_argument,       NULL, 's' },
            { "version", no_argument,       NULL, 'v' },
            { "help",    no_argument,       NULL, 'h' },
            {0, 0, 0, 0}
        };

        /* Parse the option */
        g = getopt_long(argc, argv, "q:m:svh", long_options, &option_index);

        /* We are at the end of the options */
        if (g == -1)
            break;

        /* Assign the option to variables */
        switch(g)
        {
            case 'q':
                c->query = strdup(optarg);
                break;
            case 'm':
                c->minlen = atof(optarg);
                break;
            case 's':
                c->set_match = 1;
                break;
            case 'v':
                print_version();
                return -1;
            case 'h':
                print_pileup_usage(NULL);
                return -1;
            case '?':
                sprintf(msg, "pbwtmaster [ERROR]: unknown option \"-%c\".\n", optopt);
                print_pileup_usage(msg);
                return -1;
            default:
                print_pileup_usage(NULL);
                return -1;
        }
    }

    /* Parse non-optioned arguments */
    if (optind != argc - 1)
    {
        print_pileup_usage("pbwtmaster [ERROR]: need input file name as mandatory argument");
        return -1;
    }
    else
    {
        c->instub = strdup(argv[optind]);
    }

    /* Check that a query sequence has been specified */
    if (c->query == NULL)
    {
        print_pileup_usage("pbwtmaster [ERROR]: --query option is mandatory");
        return -1;
    }

    return 0;
}

int parse_summary(int argc, char *argv[], cmd_t *c)
{
    int g = 0;
    char msg[100];

    while (1)
    {
        int option_index = 0;

        /* Declare the option table */
        static struct option long_options[] =
        {
            { "regcount", no_argument,      NULL, 'r' },
            { "version", no_argument,       NULL, 'v' },
            { "help",    no_argument,       NULL, 'h' },
            {0, 0, 0, 0}
        };

        /* Parse the option */
        g = getopt_long(argc, argv, "rvh", long_options, &option_index);

        /* We are at the end of the options */
        if (g == -1)
            break;

        /* Assign the option to variables */
        switch(g)
        {
            case 'r':
                c->reg_count = 1;
                break;
            case 'v':
                print_version();
                return -1;
            case 'h':
                print_summary_usage(NULL);
                return -1;
            case '?':
                sprintf(msg, "pbwtmaster [ERROR]: unknown option \"-%c\".\n", optopt);
                print_summary_usage(msg);
                return -1;
            default:
                print_summary_usage(NULL);
                return -1;
        }
    }

    /* Parse non-optioned arguments */
    if (optind != argc - 1)
    {
        print_summary_usage("pbwtmaster [ERROR]: need input stub as mandatory argument");
        return -1;
    }
    else
    {
        c->instub = strdup(argv[optind]);
    }

    return 0;
}

int parse_view(int argc, char *argv[], cmd_t *c)
{
    int g = 0;
    char msg[100];

    while (1)
    {
        int option_index = 0;

        /* Declare option table */
        static struct option long_options[] =
        {
            { "sites",   no_argument,       NULL, 's' },
            { "nohaps",  no_argument,       NULL, 'n' },
            { "version", no_argument,       NULL, 'v' },
            { "help",    no_argument,       NULL, 'h' },
            {0, 0, 0, 0}
        };

        /* Parse options */
        g = getopt_long(argc, argv, "snhv", long_options, &option_index);

        /* We are at the end of the options */
        if (g == -1)
        {
            break;
        }

        /* Assign option to variable */
        switch(g)
        {
            case 'n':
                c->nohaps = 1;
                break;
            case 's':
                c->only_sites = 1;
                break;
            case 'v':
                print_version();
                return -1;
            case 'h':
                print_view_usage(NULL);
                return -1;
            case '?':
                sprintf(msg, "pbwtmaster [ERROR]: unknown option \"-%c\".\n", optopt);
                print_view_usage(msg);
                return -1;
            default:
                print_view_usage(NULL);
                return -1;
        }
    }

    /* Parse non-optioned arguments */
    if (optind != argc - 1)
    {
        print_view_usage("pbwtmaster [ERROR]: need input stub as mandatory argument");
        return -1;
    }
    else
    {
        c->instub = strdup(argv[optind]);
    }

    return 0;
}

int print_main_usage(const char *msg)
{
    puts("Usage: pbwtmaster [COMMAND] [OPTION]... [INPUT FILE]\n");
    puts("Utilities for working with the PBWT format\n");
    putchar('\n');
    if (msg)
    {
        printf ("%s\n\n", msg);
    }
    puts("Commands:");
    puts("  coancesty           Construct coancestry matrix between individuals");
    puts("  convert             Convert PLINK or VCF to PBWT or vice versa");
    puts("  match               Run region matching algorithm");
    puts("  pileup              Calculate match pileup depth across chromosomes");
    puts("  summary             Produce summary of PBWT file");
    puts("  view                Dump .pbwt file to stdout");
    putchar('\n');
    return 0;
}

int print_coancestry_usage(const char *msg)
{
    puts("Usage: pbwtmaster coancestry [OPTION]... [PBWT FILE]\n");
    puts("Produce coancestry matrix for all samples in PBWT\n");
    putchar('\n');
    if (msg)
    {
        printf("%s\n\n", msg);
    }
    puts("Options:");
    puts("  --adjlist          Output graph-based adjacency list [ Default: False ]");
    puts("  --diploid          Output diploid rather than haploid-based measures");
    puts("  --set              Find only set-maximal matches [ Default: all matches ]");
    puts("  --sites            Print site indices [ Default: false ]");
    puts("  --count            Coancestry matrix will have match count [ Default: total length ]");
    puts("  --minlen   FLOAT   Minimum match size (cM) [ Default: 0.5 cM ]");
    puts("  --version          Print version number and exit");
    puts("  --help             Display this help message and exit");
    putchar('\n');
    return 0;
}

int print_convert_usage(const char *msg)
{
    puts("Usage: pbwtmaster convert [OPTION]... [INPUT STUB]\n");
    puts("Convert PLINK or VCF to PBWT or vice versa\n");
    putchar('\n');
    if (msg)
    {
        printf("%s\n\n", msg);
    }
    puts("Options:");
    puts("  --phased            Input data are phased");
    puts("  --vcf               Input file is VCF format");
    puts("  --map     <FILE>    Popmap file (for VCF input)");
    puts("  --reg               Input PLINK stub includes a REG file");
    puts("  --out      <STR>    Output stub (.pbwt extension will be added)");
    puts("  --query    <STR>    Use only this region/population (requires -r switch)");
    puts("  --version           Print version number and exit");
    puts("  --help              Display this help message and exit");
    putchar('\n');
    return 0;
}

int print_match_usage(const char *msg)
{
    puts("Usage: pbwtmaster match [OPTION]... [PBWT FILE]\n");
    puts("Retrieve match sets in PBWT using query\n");
    putchar('\n');
    if (msg)
    {
        printf("%s\n\n", msg);
    }
    puts("Options:");
    puts("  --minlen   FLOAT   Minimum match size (cM) [ Default: 0.5 cM ]");
    puts("  --query    STR     String identifier of haplotypes to mark as query");
    puts("  --all              Print a list of all individual matches with query");
    puts("  --set              Find only set-maximal matches [ Default: all matches ]");
    puts("  --sites            Print site indices [ Default: false ]");
    puts("  --version          Print version number and exit");
    puts("  --help             Display this help message and exit");
    putchar('\n');
    return 0;
}

int print_pileup_usage(const char *msg)
{
    puts("Usage: pbwtmaster pileup [OPTION]... [PBWT FILE]\n");
    puts("Calculate match pileup depth across chromosomes\n");
    putchar('\n');
    if (msg)
    {
        printf("%s\n\n", msg);
    }
    puts("Options:");
    puts("  --minlen   FLOAT   Minimum match size (cM) [ Default: 0.5 cM ]");
    puts("  --query    STR     String identifier of haplotypes to mark as query");
    puts("  --set              Find only set-maximal matches [ Default: all matches ]");
    puts("  --version          Print version number and exit");
    puts("  --help             Display this help message and exit");
    putchar('\n');
    return 0;
}

int print_summary_usage(const char *msg)
{
    puts("Usage: pbwtmaster summary [OPTION]... [INPUT STUB]\n");
    puts("Print summary describing .pbwt file\n");
    putchar('\n');
    if (msg)
    {
        printf ("%s\n\n", msg);
    }
    puts("Options:");
    puts("  --regcount          Print region sizes");
    puts("  --version           Print version number and exit");
    puts("  --help              Display this help message and exit");
    putchar('\n');
    return 0;
}

int print_view_usage(const char *msg)
{
    puts("Usage: pbwtmaster view [OPTION]... [PBWT FILE]\n");
    puts("Print contents of .pbwt file to stdout\n");
    putchar('\n');
    if (msg)
    {
        printf("%s\n\n", msg);
    }
    puts("Options:");
    puts("  --nohaps            Omit haplotype states-- only print sample metadata");
    puts("  --sites             Print only site information");
    puts("  --version           Print version number and exit");
    puts("  --help              Display this help message and exit");
    putchar('\n');
    return 0;
}

void print_version(void)
{
    printf("pbwtmaster: %s\n", Version);
    printf("libpbwt: %s\n", libpbwt_version());
    printf("libplink_lite: %s\n", plink_version());
    printf("htslib: %s\n", hts_version());
}
