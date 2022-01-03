#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pbwtutil.h"

void report_adjlist(pbwt_t *b, const size_t first, const size_t second, const size_t begin, const size_t end)
{
    printf("%s\t%s\t%1.4lf\t%s\t%s\n", b->sid[first], b->sid[second],
           b->cm[end] - b->cm[begin], b->reg[first], b->reg[second]);
}

void report_adjlist_with_sites(pbwt_t *b, const size_t first, const size_t second, const size_t begin, const size_t end)
{
    printf("%s\t%s\t%1.4lf\t%s\t%s\t%zu\t%zu\n", b->sid[first], b->sid[second],
            b->cm[end] - b->cm[begin], b->reg[first], b->reg[second], begin, end);
}

void add_interval(pbwt_t *b, const size_t first, const size_t second, const size_t begin, const size_t end)
{
    b->intree = match_insert(b->intree, first, second, begin, end);
}


void add_nmatch(pbwt_t *b, const size_t first, const size_t second, const size_t begin, const size_t end)
{
	if (b->nmatrix == NULL)
	{
		b->nmatrix = (size_t **)malloc(b->nsam * sizeof(size_t *));
		if (b->nmatrix == NULL)
		{
			return;
		}
		size_t i = 0;
		for (i = 0; i < b->nsam; ++i)
		{
			b->nmatrix[i] = (size_t *)calloc(b->nsam, sizeof(size_t));
			if (b->nmatrix[i] == NULL)
			{
				return;
			}
		}
	}
	b->nmatrix[first][second]++;
	b->nmatrix[second][first]++;
}

void add_coancestry(pbwt_t *b, const size_t first, const size_t second, const size_t begin, const size_t end)
{
	if (b->cmatrix == NULL)
	{
		b->cmatrix = (double **)malloc(b->nsam * sizeof(double*));
		if (b->cmatrix == NULL)
		{
			return;
		}
		size_t i = 0;
		for (i = 0; i < b->nsam; i++)
		{
			b->cmatrix[i] = (double *)calloc(b->nsam, sizeof(double));
			if (b->cmatrix[i] == NULL)
			{
				return;
			}
		}
	}
	double length = b->cm[end] - b->cm[begin];
	b->cmatrix[first][second] += length;
	b->cmatrix[second][first] = b->cmatrix[first][second];
}

void add_region(pbwt_t *b, const size_t first, const size_t second, const size_t begin, const size_t end)
{
	if (b->reghash == NULL)
	{
		b->reghash = kh_init(floats);
	}
	int a = 0;
	size_t qs = 0;
	khint_t k = 0;
	double length = b->cm[end] - b->cm[begin];
	qs = b->is_query[first] ? second : first;
	k = kh_put(floats, b->reghash, b->reg[qs], &a);
	if (a == 0)
	{
		double ent = kh_value(b->reghash, k);
		ent += length;
		kh_value(b->reghash, k) = ent;
	}
	else
	{
		kh_value(b->reghash, k) = length;
	}
}
