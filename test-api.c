#include <stdio.h>
#include <stdlib.h>
#include "plink.h"

#define FALSE 0
#define TRUE 1

int main (int argc, char *argv[])
{
    /*size_t i = 0;*/
    plink_t *p = read_plink(argv[1], TRUE);
    if (p == NULL)
        return 1;

    /* How to use query_reg() function */
    char q[] = "1651_Sir43";
    char *res = query_reg(q, p->reg);
    fprintf (stderr, "Query: %s\nResult:%s\n", q, res);

    return 0;
}