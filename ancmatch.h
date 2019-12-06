#ifndef ANCMATCH_H
#define ANCMATCH_H

#include <htslib/khash.h>
#include <htslib/vcf.h>
#include <pbwt.h>

/* Define mode mappings */

enum Mode {COANCESTRY, CONVERT, REPORT, RUN, VIEW};


/* Define data structures */

typedef struct cmdl
{
    enum Mode mode;
    int with_vcf;
    int has_reg;
    int is_phased;
    double minlen;
    char *popmap;
    char *outfile;
    char *query;
    char *instub;
    int (*mode_func)(struct cmdl *);
} cmd_t;


/* Function prototypes */

extern cmd_t *parse_cmdl(int argc, char *argv[]);

extern int pbwt_coancestry(cmd_t *);

extern int pbwt_convert_plink(cmd_t *);

extern int pbwt_convert_vcf(cmd_t *);

extern int pbwt_report(cmd_t *);

extern int pbwt_run(cmd_t *);

extern int pbwt_view(cmd_t *);

#endif

