#include <stdio.h>
#include <stdlib.h>
#include "ancmatch.h"


int pbwt_view(cmd_t *c)
{
    int v = 0;
    pbwt_t *b = NULL;

    b = pbwt_read(c->instub);
    v = pbwt_uncompress(b);
    v = pbwt_print(b);
    pbwt_destroy(b);

    return v;
}
