#include "plink.h"

bed_t *
bed_read (const char *bedfile, uint64_t nsam, uint64_t nsnp, unsigned char *data)
{
    bed_t *bed;
    FILE *fin;

    /* Allocate heap memory for bed_t structure */
    bed = (bed_t *) malloc (sizeof(bed_t));
    if (bed == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Open binary input stream */
    fin = fopen (bedfile, "rb");
    if (fin == NULL)
    {
        perror ("libplink [ERROR]");
        return NULL;
    }

    /* Read first three bytes */
    const int h1 = fgetc (fin);
    const int h2 = fgetc (fin);
    const int smo = fgetc (fin);

    /* No data */
    if (feof (fin) != 0)
    {
        fputs ("libplink [ERROR]: input bed/hap appears to be truncated\n", stderr);
        return NULL;
    }

    /* Determine memory block size based on orientation */
    const size_t record_size = smo ? bed_record_size(nsam) : bed_record_size(nsnp);
    const size_t expected_size = smo ? (nsnp * record_size) : (nsam * record_size);

    /* Allocate memory for data blob */
    if (data == NULL)
    {
        data = (unsigned char *) malloc (expected_size * sizeof(unsigned char));
        if (data == NULL)
            return NULL;
    }

    /* Read data blob into memory */
    const size_t bytesread = fread (data, sizeof(unsigned char), expected_size, fin);

    /* Possibly truncated file */
    if (bytesread != expected_size)
    {
        fputs ("Unexpected EOF reading binary data", stderr);
        free (data);
        return NULL;
    }

    /* Check if we are at the end of the file */
    fgetc (fin);
    if (feof (fin) == 0)
    {
        fputs ("Binary ped file larger than expected", stderr);
        free (data);
        return NULL;
    }

    /* Populate bed data structure */
    bed->header1 = h1;
    bed->header2 = h2;
    bed->orientation = smo;
    bed->record_size = record_size;
    bed->size = bytesread;
    bed->data = data;
    bed->phased = (bed->header2 == HAP_MAGIC2) ? TRUE : FALSE;

    return bed;
}

void
bed_destroy (bed_t *bed)
{
    if (bed != NULL)
    {
        if (bed->data != NULL)
            free (bed->data);
        free (bed);
    }

    return;
}
