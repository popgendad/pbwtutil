#include "ancmatch.h"

int main(int argc, char **argv)
{
    cmd_t *c = parse_cmdl(argc, argv);
    if (!c)
    {
        return -1;
    }
    else
    {
        return (*c->mode_func)(c);
    }
}
