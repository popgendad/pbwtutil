#include "pbwtmaster.h"

int main(int argc, char **argv)
{
    cmd_t *c = NULL;
 
    c = parse_args(argc, argv);
    if (!c)
    {
        return -1;
    }

    return (*c->mode_func)(c);
}
