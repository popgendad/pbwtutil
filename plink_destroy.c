#include "plink.h"

void
plink_destroy (plink_t *p)
{
    if (p->bed != NULL)
        bed_destroy (p->bed);
    if (p->bim != NULL)
        bim_destroy (p->bim, p->nsnp);
    if (p->fam != NULL)
        fam_destroy (p->fam, p->nsam);
    if (p->reg != NULL)
        reg_destroy (p->reg, p->nsam);
    if (p->bim_index != NULL)
        kh_destroy(integer, p->bim_index);
    if (p->fam_index != NULL)
        kh_destroy(integer, p->fam_index);
    if (p->reg_index != NULL)
        kh_destroy(integer, p->reg_index);
    free (p);

    return;
}
