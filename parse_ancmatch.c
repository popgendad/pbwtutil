#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "ancmatch.h"

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
int parse_report(int, char **, cmd_t *);
int parse_run(int, char **, cmd_t *);
int parse_view(int, char **, cmd_t *);
int print_main_usage(const char *);
int print_coancestry_usage(const char *);
int print_convert_usage(const char *);
int print_report_usage(const char *);
int print_view_usage(const char *);
int print_run_usage(const char *);

cmd_t *parse_cmdl(int argc, char *argv[])
{
    char *mode = NULL;
    cmd_t *c = NULL;
    int (*parse_func)(int, char **, cmd_t *);

    opterr = 0;

    /* Allocate memory for command line data structure */
    c = (cmd_t *)malloc(sizeof(cmd_t));
    if (c == NULL)
    {
        fputs("ancmatch [ERROR]: memory allocation failure", stderr);
        return NULL;
    }

    /* Initialize default values for command line options */
    c->has_reg = 0;
    c->is_phased = 0;
    c->with_vcf = 0;
    c->nohaps = 0;
    c->minlen = 0.5;

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
    else if (strcmp(mode, "report") == 0)
    {
        c->mode = REPORT;
        c->mode_func = &pbwt_report;
        parse_func = &parse_report;
    }
    else if (strcmp(mode, "run") == 0)
    {
        c->mode = RUN;
        c->mode_func = &pbwt_run;
        parse_func = &parse_run;
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
            sprintf(msg, "ancmatch [ERROR]: unknown command \'%s\'", mode);
            print_main_usage(msg);
        }
        else
        {
            print_main_usage("ancmatch [ERROR]: need to specify a command");
        }

        return NULL;
    }

    /* Call the command-specific parsing function */
    int r = (*parse_func)(argc, argv, c);
    if (r)
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
            { "minlen",  required_argument, NULL, 'm' },
            { "version", no_argument,       NULL, 'v' },
            { "help",    no_argument,       NULL, 'h' },
            {0, 0, 0, 0}
        };

        /* Parse the option */
        g = getopt_long(argc, argv, "vhm:", long_options, &option_index);

        /* We are at the end of the options */
        if (g == -1)
        {
            break;
        }

        /* Assign the option to variables */
        switch(g)
        {
            case 'm':
                c->minlen = atof(optarg);
                break;
            case 'v':
                printf("ancmatch: %s\n", Version);
                printf("libpbwt: %s\n", pbwt_version());
                return 1;
            case 'h':
                print_coancestry_usage(NULL);
                return 1;
            case '?':
                sprintf(msg, "ancmatch [ERROR]: unknown option \"-%c\".\n", optopt);
                print_coancestry_usage(msg);
                return 1;
            default:
                print_coancestry_usage(NULL);
                return 1;
        }
    }

    /* Parse non-optioned arguments */
    if (optind != argc - 1)
    {
        print_coancestry_usage("ancmatch [ERROR]: need input file name as mandatory argument");
        return 1;
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
                printf("ancmatch: %s\n", Version);
                printf("libpbwt: %s\n", pbwt_version());
                return 1;
            case 'h':
                print_convert_usage(NULL);
                return 1;
            case '?':
                sprintf(msg, "ancmatch [ERROR]: unknown option \"-%c\".\n", optopt);
                print_convert_usage(msg);
                return 1;
            default:
                print_convert_usage(NULL);
                return 1;
        }
    }

    /* Parse non-optioned arguments */
    if (optind != argc - 1)
    {
        print_convert_usage("ancmatch [ERROR]: need input stub as mandatory argument");
        return 1;
    }
    else
    {
        c->instub = strdup(argv[optind]);
    }

    /* --out switch is actually mandatory */
    if (!c->outfile)
    {
        print_convert_usage("ancmatch [ERROR]: --out <STR> is a mandatory argument for convert");
        return 1;
    }

    return 0;
}

int parse_report(int argc, char *argv[], cmd_t *c)
{
    int g = 0;
    char msg[100];

    while (1)
    {
        int option_index = 0;

        /* Declare the option table */
        static struct option long_options[] =
        {
            { "version", no_argument,       NULL, 'v' },
            { "help",    no_argument,       NULL, 'h' },
            {0, 0, 0, 0}
        };

        /* Parse the option */
        g = getopt_long (argc, argv, "vh", long_options, &option_index);

        /* We are at the end of the options */
        if (g == -1)
            break;

        /* Assign the option to variables */
        switch(g)
        {
            case 'v':
                printf("ancmatch: %s\n", Version);
                printf("libpbwt: %s\n", pbwt_version());
                return 1;
            case 'h':
                print_report_usage(NULL);
                return 1;
            case '?':
                sprintf(msg, "ancmatch [ERROR]: unknown option \"-%c\".\n", optopt);
                print_report_usage(msg);
                return 1;
            default:
                print_report_usage(NULL);
                return 1;
        }
    }

    /* Parse non-optioned arguments */
    if (optind != argc - 1)
    {
        print_report_usage ("ancmatch [ERROR]: need input stub as mandatory argument");
        return 1;
    }
    else
    {
        c->instub = strdup(argv[optind]);
    }

    return 0;
}

int parse_run(int argc, char *argv[], cmd_t *c)
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
            { "version", no_argument,       NULL, 'v' },
            { "help",    no_argument,       NULL, 'h' },
            {0, 0, 0, 0}
        };

        /* Parse the option */
        g = getopt_long(argc, argv, "vhq:m:", long_options, &option_index);

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
            case 'm':
                c->minlen = atof(optarg);
                break;
            case 'v':
                printf("ancmatch: %s\n", Version);
                printf("libpbwt: %s\n", pbwt_version());
                return 1;
            case 'h':
                print_run_usage(NULL);
                return 1;
            case '?':
                sprintf(msg, "ancmatch [ERROR]: unknown option \"-%c\".\n", optopt);
                print_run_usage(msg);
                return 1;
            default:
                print_run_usage(NULL);
                return 1;
        }
    }

    /* Parse non-optioned arguments */
    if (optind != argc - 1)
    {
        print_run_usage("ancmatch [ERROR]: need input file name as mandatory argument");
        return 1;
    }
    else
    {
        c->instub = strdup(argv[optind]);
    }

    /* Check that a query sequence has been specified */
    if (c->query == NULL)
    {
        print_run_usage("ancmatch [ERROR]: --query option is mandatory");
        return 1;
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
            { "nohaps",  no_argument,       NULL, 'n' },
            { "version", no_argument,       NULL, 'v' },
            { "help",    no_argument,       NULL, 'h' },
            {0, 0, 0, 0}
        };

        /* Parse options */
        g = getopt_long(argc, argv, "hv", long_options, &option_index);

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
            case 'v':
                printf("ancmatch: %s\n", Version);
                printf("libpbwt: %s\n", pbwt_version());
                return 0;
            case 'h':
                print_view_usage(NULL);
                return 1;
            case '?':
                sprintf(msg, "ancmatch [ERROR]: unknown option \"-%c\".\n", optopt);
                print_view_usage(msg);
                return 1;
            default:
                print_view_usage(NULL);
                return 1;
        }
    }

    /* Parse non-optioned arguments */
    if (optind != argc - 1)
    {
        print_view_usage("ancmatch [ERROR]: need input stub as mandatory argument");
        return 1;
    }
    else
    {
        c->instub = strdup(argv[optind]);
    }

    return 0;
}

int print_main_usage(const char *msg)
{
    puts("Usage: ancmatch [COMMAND] [OPTION]... [INPUT FILE]\n");
    puts("Utilities for working with the PBWT format\n");
    putchar('\n');
    if (msg)
        printf ("%s\n\n", msg);
    puts("Commands:");
    puts("  coancesty           Construct coancestry matrix between individuals");
    puts("  convert             Convert PLINK or VCF to PBWT or vice versa");
    puts("  report              Report on compression available in pbwt");
    puts("  run                 Run region matching algorithm");
    puts("  view                Dump .pbwt file to stdout");
    putchar('\n');
    return 0;
}

int print_coancestry_usage(const char *msg)
{
    puts("Usage: ancmatch coancestry [OPTION]... [PBWT FILE]\n");
    puts("Produce coancestry matrix for all samples in PBWT\n");
    putchar('\n');
    if (msg)
        printf("%s\n\n", msg);
    puts("Options:");
    puts("  --minlen   FLOAT   Minimum match size (cM) [ Default: 0.5 cM ]");
    puts("  --version          Print version number and exit");
    puts("  --help             Display this help message and exit");
    putchar('\n');
    return 0;
}

int print_convert_usage(const char *msg)
{
    puts("Usage: ancmatch convert [OPTION]... [INPUT STUB]\n");
    puts("Convert PLINK or VCF to PBWT or vice versa\n");
    putchar('\n');
    if (msg)
        printf("%s\n\n", msg);
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

int print_report_usage(const char *msg)
{
    puts("Usage: ancmatch report [OPTION]... [INPUT STUB]\n");
    puts("Print report describing pbwt file\n");
    putchar('\n');
    if (msg)
        printf ("%s\n\n", msg);
    puts("Options:");
    puts("  --version           Print version number and exit");
    puts("  --help              Display this help message and exit");
    putchar('\n');
    return 0;
}

int print_run_usage(const char *msg)
{
    puts("Usage: ancmatch run [OPTION]... [PBWT FILE]\n");
    puts("Retrieve set-maximal matches in PBWT using query\n");
    putchar('\n');
    if (msg)
        printf("%s\n\n", msg);
    puts("Options:");
    puts("  --minlen   FLOAT   Minimum match size (cM) [ Default: 0.5 cM ]");
    puts("  --query    STR     String identifier of haplotypes to mark as query");
    puts("  --version          Print version number and exit");
    puts("  --help             Display this help message and exit");
    putchar('\n');
    return 0;
}

int print_view_usage(const char *msg)
{
    puts("Usage: ancmatch view [OPTION]... [PBWT FILE]\n");
    puts("Print contents of .pbwt file to stdout\n");
    putchar('\n');
    if (msg)
        printf("%s\n\n", msg);
    puts("Options:");
    puts("  --nohaps            Omit haplotype states-- only print sample metadata");
    puts("  --version           Print version number and exit");
    puts("  --help              Display this help message and exit");
    putchar('\n');
    return 0;
}

