#include "plink.h"

fam_t *
fam_read (const char *famfile, size_t *nl)
{
    size_t lines_read;
    size_t total_alloc = CHUNK_SIZE;
    char line[LINE_LENGTH];
    const char delim[] = " \t\n";
    fam_t *fam;
    FILE *fin;

    /* Open fam input file stream */
    fin = fopen (famfile, "r");
    if (fin == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Allocate heap memory for fam_t structure */
    fam = (fam_t *) malloc (total_alloc * sizeof(fam_t));
    if (fam == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Initialize line counter */
    lines_read = 0;

    /* Iterate through the lines of the fam file */
    while (fgets (line, LINE_LENGTH, fin) != NULL)
    {
        ++lines_read;
        if (lines_read > total_alloc)
        {
            total_alloc += CHUNK_SIZE;
            fam = (fam_t *) realloc (fam, total_alloc * sizeof(fam_t));
        }

        /* Parse the fam entry */
        char *token;
        size_t i;

        /* Define index of position */
        i = lines_read - 1;

        /* Parse the tokens */
        token  = strtok (line, delim);
        fam[i].fid = strdup (token);
        token = strtok (NULL, delim);
        fam[i].iid = strdup (token);
        token = strtok (NULL, delim);
        fam[i].pid = strdup (token);
        token = strtok (NULL, delim);
        fam[i].mid = strdup (token);
        token = strtok (NULL, delim);
        fam[i].sex = strdup (token);
        token = strtok (NULL, delim);
        fam[i].phe = strdup (token);
    }

    /* Re-size the array holding fam data */
    fam = (fam_t *) realloc (fam, lines_read * sizeof(fam_t));

    /* Close the input file stream */
    fclose (fin);

    /* Assign the number of lines read from the fam file */
    *nl = lines_read;

    return fam;
}


int
fam_write (const char *outfile, const fam_t *fam, const size_t nl)
{
    size_t i;
    FILE *fout;

    fout = fopen (outfile, "w");
    if (fout == NULL)
    {
        perror ("libplink [ERROR]");
        return -1;
    }

    for (i = 0; i < nl; ++i)
        fprintf (fout, "%s\t%s\t%s\t%s\t%s\t%s\n",
                 fam[i].fid, fam[i].iid, fam[i].pid, fam[i].mid,
                 fam[i].sex, fam[i].phe);

    fclose (fout);

    return 0;
}


khash_t(integer) *
fam_index (const fam_t *fam, const size_t nl)
{
    int a;
    size_t i;
    khint_t k;
    khash_t(integer) *fx;

    /* Initialize fam hash */
    fx = kh_init(integer);

    /* Iterate through lines in bim and store indices keyed on sample identifier */
    for (i = 0; i < nl; ++i)
    {
        k = kh_put(integer, fx, fam[i].iid, &a);
        kh_value(fx, k) = i;
    }

    return fx;
}


void
fam_destroy (fam_t *fam, const size_t nsam)
{
    size_t i;

    if (fam != NULL)
    {
        for (i = 0; i < nsam; ++i)
        {
            if (fam[i].fid != NULL)
                free (fam[i].fid);
            if (fam[i].iid != NULL)
                free (fam[i].iid);
            if (fam[i].pid != NULL)
                free (fam[i].pid);
            if (fam[i].mid != NULL)
                free (fam[i].mid);
            if (fam[i].sex != NULL)
                free (fam[i].sex);
            if (fam[i].phe != NULL)
                free (fam[i].phe);
        }
        free (fam);
    }

    return;
}
