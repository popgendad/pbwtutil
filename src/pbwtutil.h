#ifndef PBWTUTIL_H
#define PBWTUTIL_H

#include <htslib/khash.h>
#include <htslib/vcf.h>
#include <pbwt.h>

/* Define mode mappings */

enum Mode {COANCESTRY, CONVERT, MATCH, PILEUP, SUMMARY, VIEW};


/* Define data structures */

typedef struct cmdl
{
    enum Mode mode;
    int with_vcf;
    int has_reg;
    int is_phased;
    int nohaps;
    int match_all;
    int print_sites;
    int only_sites;
    int count_only;
    int reg_count;
    int adjlist;
    int out_diploid;
    int set_match;
    double minlen;
    char *popmap;
    char *outfile;
    char *query;
    char *instub;
    int (*mode_func)(const struct cmdl *);
} cmd_t;


/* Function prototypes */

extern cmd_t *parse_args(int argc, char *argv[]);

extern int pbwt_coancestry(const cmd_t *);

extern int pbwt_convert_plink(const cmd_t *);

extern int pbwt_convert_vcf(const cmd_t *);

extern int pbwt_match(const cmd_t *);

extern int pbwt_pileup(const cmd_t *);

extern int pbwt_summary(const cmd_t *);

extern int pbwt_view(const cmd_t *);

extern void add_interval(pbwt_t *, const size_t, const size_t, const size_t, const size_t);

extern  void report_adjlist(pbwt_t *, const size_t, const size_t, const size_t, const size_t);

extern void report_adjlist_with_sites(pbwt_t *, const size_t, const size_t, const size_t, const size_t);

extern void add_coancestry(pbwt_t *, const size_t, const size_t, const size_t, const size_t);

extern void add_nmatch(pbwt_t *, const size_t, const size_t, const size_t, const size_t);

extern void add_region(pbwt_t *, const size_t, const size_t, const size_t, const size_t);

#endif
