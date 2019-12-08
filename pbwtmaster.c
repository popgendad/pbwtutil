#include "pbwtmaster.h"

int main(int argc, char **argv)
{
    cmd_t *c = parse_args(argc, argv);
    if (!c)
    {
        return -1;
    }
    else
    {
        return (*c->mode_func)(c);
    }
}
