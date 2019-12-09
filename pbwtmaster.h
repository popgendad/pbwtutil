#ifndef PBWTMASTER_H
#define PBWTMASTER_H

#include <htslib/khash.h>
#include <htslib/vcf.h>
#include <pbwt.h>

/* Define mode mappings */

enum Mode {COANCESTRY, CONVERT, MATCH, SUMMARY, VIEW};


/* Define data structures */

typedef struct cmdl
{
    enum Mode mode;
    int with_vcf;
    int has_reg;
    int is_phased;
    int nohaps;
    int match_all;
    double minlen;
    char *popmap;
    char *outfile;
    char *query;
    char *instub;
    int (*mode_func)(struct cmdl *);
} cmd_t;


/* Function prototypes */

extern cmd_t *parse_args(int argc, char *argv[]);

extern int pbwt_coancestry(cmd_t *);

extern int pbwt_convert_plink(cmd_t *);

extern int pbwt_convert_vcf(cmd_t *);

extern int pbwt_match(cmd_t *);

extern int pbwt_summary(cmd_t *);

extern int pbwt_view(cmd_t *);

#endif

