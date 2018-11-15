#include "plink.h"


reg_t *
reg_read (const char *regfile, size_t *nl)
{
    size_t lines_read;
    size_t total_alloc = CHUNK_SIZE;
    reg_t *reg;
    FILE *fin;
    char line[LINE_LENGTH];
    const char delim[] = " \t\n";

    /* Open the input reg file stream */
    fin = fopen (regfile, "r");
    if (fin == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    reg = (reg_t *) malloc (total_alloc * sizeof(reg_t));
    if (reg == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Initialize the line counter */
    lines_read = 0;

    /* Iterate through the lines of the reg file */
    while (fgets (line, LINE_LENGTH, fin) != NULL)
    {
        ++lines_read;
        if (lines_read > total_alloc)
        {
            total_alloc += CHUNK_SIZE;
            reg = (reg_t *) realloc (reg, total_alloc * sizeof(reg_t));
        }

        /* Parse the reg entry */
        char *token;
        size_t i;

        /* Define the index of position */
        i = lines_read - 1;

        /* Parse the tokens */
        token = strtok (line, delim);
        reg[i].fid = strdup (token);
        token = strtok (NULL, delim);
        reg[i].iid = strdup (token);
        token = strtok (NULL, delim);
        reg[i].pid = strdup (token);
        token = strtok (NULL, delim);
        reg[i].mid = strdup (token);
        token = strtok (NULL, delim);
        reg[i].sex = strdup (token);
        token = strtok (NULL, delim);
        reg[i].phe = strdup (token);
        token = strtok (NULL, delim);
        reg[i].pop = strdup (token);
        token = strtok (NULL, delim);
        reg[i].reg = strdup (token);
    }

    /* Re-size the heap memory for the reg structure */
    reg = (reg_t *) realloc (reg, lines_read * sizeof(reg_t));

    /* Close the input file stream */
    fclose (fin);

    /* Assign the number of lines read from the reg file */
    *nl = lines_read;

    return reg;
}


int
reg_write (const char *outfile, const reg_t *reg, const size_t nl)
{
    size_t i;
    FILE *fout;

    /* Open output file stream */
    fout = fopen (outfile, "w");
    if (fout == NULL)
    {
        perror ("libplink [ERROR]");
        return -1;
    }

    /* Print entry to file */
    for (i = 0; i < nl; ++i)
        fprintf (fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                 reg[i].fid, reg[i].iid, reg[i].pid, reg[i].mid,
                 reg[i].sex, reg[i].phe, reg[i].pop, reg[i].reg);

    /* Close the output file stream */
    fclose (fout);

    return 0;
}


khash_t(integer) *
reg_index (const reg_t *reg, const size_t nl)
{
    int a;
    size_t i;
    khint_t k;
    khash_t(integer) *rx;

    /* Initialize reg file hash */
    rx = kh_init(integer);

    /* Load reg data into hash */
    for (i = 0; i < nl; ++i)
    {
        k = kh_put(integer, rx, reg[i].iid, &a);
        kh_value(rx, k) = i;
    }

    return rx;
}


void
reg_destroy (reg_t *reg, const size_t nsam)
{
    size_t i;

    if (reg != NULL)
    {
        for (i = 0; i < nsam; ++i)
        {
            if (reg[i].fid != NULL)
                free (reg[i].fid);
            if (reg[i].iid != NULL)
                free (reg[i].iid);
            if (reg[i].pid != NULL)
                free (reg[i].pid);
            if (reg[i].mid != NULL)
                free (reg[i].mid);
            if (reg[i].sex != NULL)
                free (reg[i].sex);
            if (reg[i].phe != NULL)
                free (reg[i].phe);
            if (reg[i].pop != NULL)
                free (reg[i].pop);
            if (reg[i].reg != NULL)
                free (reg[i].reg);
        }
        free (reg);
    }

    return;
}
