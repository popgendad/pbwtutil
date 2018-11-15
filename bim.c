#include "plink.h"

/* Unknown allele value in .bim file */
const char UNKNOWN_VARIANT = '0';


bim_t *
bim_read (const char *bimfile, size_t *nl)
{
    size_t lines_read;
    size_t total_alloc = CHUNK_SIZE;
    char line[LINE_LENGTH];
    const char delim[] = " \t\n";
    bim_t *bim;
    FILE *fin;

    /* Open input file stream */
    fin = fopen (bimfile, "r");
    if (fin == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Initialize line counter */
    lines_read = 0;

    /* Allocate heap memory for bim_t structure */
    bim = (bim_t *) malloc (total_alloc * sizeof(bim_t));
    if (bim == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Iterate through lines of bim file */
    while (fgets (line, LINE_LENGTH, fin) != NULL)
    {
        ++lines_read;
        if (lines_read > total_alloc)
        {
            total_alloc += CHUNK_SIZE;
            bim = (bim_t *) realloc (bim, total_alloc * sizeof(bim_t));
        }

        /* Parse the bim entry */
        char *token;
        size_t i;

        /* Define index position */
        i = lines_read - 1;
   
        /* Parse tokens */
        token = strtok (line, delim);
        if (strcmp (token, "X") == 0)
            bim[i].ch = 23;
        else if (strcmp (token, "Y") == 0)
            bim[i].ch = 24;
        else if (strcmp (token, "XY") == 0)
            bim[i].ch = 25;
        else if (strcmp (token, "MT") == 0)
            bim[i].ch = 26;
        else
            bim[i].ch = atoi (token);
        token = strtok (NULL, delim);
        bim[i].rsid = strdup (token);
        token = strtok (NULL, delim);
        bim[i].cM = strtod (token, NULL);
        token = strtok (NULL, delim);
        bim[i].bp = strtoul (token, NULL, 0);
        token = strtok (NULL, delim);
        bim[i].a0 = token[0];
        token = strtok (NULL, delim);
        bim[i].a1 = token[0];
    }

    /* Resize bim array heap memory */
    bim = (bim_t *) realloc (bim, lines_read * sizeof(bim_t));

    /* Close the input file stream */
    fclose (fin);

    /* Set total number of lines read from bim file */
    *nl = lines_read;

    return bim;
}

int
bim_write (const char *outfile, const bim_t *bim, const size_t nl)
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

    for (i = 0; i < nl; ++i)
        fprintf (fout, "%d\t%s\t%lf\t%lu\t%c\t%c\n",
                 bim[i].ch, bim[i].rsid, bim[i].cM, bim[i].bp,
                 bim[i].a0, bim[i].a1);

    /* Close output file stream */
    fclose (fout);

    return 0;
}


khash_t(integer) *
bim_index (const bim_t *bim, const size_t nl)
{
    int a;
    size_t i;
    khint_t k;
    khash_t(integer) *bx;

    /* Initialize bim hash */
    bx = kh_init(integer);

    /* Iterate through lines in bim and store indices keyed on rsid */
    for (i = 0; i < nl; ++i)
    {
        k = kh_put(integer, bx, bim[i].rsid, &a);
        kh_value(bx, k) = i;
    }

    return bx;
}


void
bim_destroy (bim_t *bim, const size_t nsnp)
{
    size_t i;

    if (bim != NULL)
    {
        for (i = 0; i < nsnp; ++i)
            if (bim[i].rsid != NULL)
                free (bim[i].rsid);
        free (bim);
    }

    return;
}
