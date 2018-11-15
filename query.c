#include "plink.h"

char *
query_reg (plink_t *p, const char *iid)
{
    uint64_t i;
    char *result;
    khint_t k;

    /* Query reg hash for sample identifier */
    k = kh_get(integer, p->reg_index, iid);

    /* Store resulting reg entry */
    if (k != kh_end(p->reg_index))
    {
        i = kh_value(p->reg_index, k);
        result = strdup (p->reg[i].reg);
        return result;
    }
    else
        return NULL;
}


char *
query_pop (plink_t *p, const char *iid)
{
    uint64_t i;
    char *result;
    khint_t k;

    /* Query reg hash for sample identifier */
    k = kh_get(integer, p->reg_index, iid);

    /* Store resulting reg entry */
    if (k != kh_end(p->reg_index))
    {
        i = kh_value(p->reg_index, k);
        result = strdup (p->reg[i].pop);
        return result;
    }
    else
        return NULL;
}
