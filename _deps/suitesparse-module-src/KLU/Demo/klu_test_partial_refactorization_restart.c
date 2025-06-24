/* klu_test_partial_refactorization_restart: a partial refactorization KLU demo, for testing */

#include <stdio.h>
#include <math.h>
#include "klu.h"

#define TOLERANCE 1e-8

int    n = 10 ;
int    Ap [ ] = { 0,  2,  3,  6,  9, 12, 15, 20, 21, 27, 31 } ;
int    Ai [ ] = { 0, 8, 1, 2, 6, 9, 3, 4, 6, 4, 5, 8, 4, 5, 8, 2, 3, 6, 8, 9, 7, 0, 4, 5, 6, 8, 9, 2, 6, 8, 9 } ;
double Ax [ ] = {8.18413247, 0.31910091, 0.95960852, 7.9683539 , 3.27076739,
       9.3203983 , 2.94765012, 0.41596915, 8.55865174, 3.26336244,
       2.56358029, 7.29705002, 9.42558416, 6.80016439, 5.82804034,
       9.39211732, 9.31241378, 0.35525264, 7.68775477, 5.48634592,
       2.80075036, 2.36812029, 1.13390547, 9.71284119, 6.02692506,
       4.03715243, 4.36857613, 0.54369597, 6.86482384, 6.46735381,
       4.76819917 } ;
double Ax_new [ ] = {8.18413247, 0.31910091, 0.95960852, 7.9683539 , 3.27076739,
       9.3203983 , 2.94765012, 1.41596915, 8.55865174, 3.26336244,
       2.56358029, 7.29705002, 9.42558416, 6.80016439, 5.82804034,
       9.39211732, 9.31241378, 0.35525264, 7.68775477, 5.48634592,
       2.80075036, 2.36812029, 1.13390547, 9.71284119, 7.02692506,
       4.03715243, 4.36857613, 0.54369597, 6.86482384, 6.46735381,
       4.76819917 } ;
double b [ ] = {1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0} ;
double c [ ] = {1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0} ;

int main (void)
{
    klu_symbolic *Symbolic ;
    klu_numeric *Numeric ;
    klu_common Common ;
    double max_error, error = 0.0;
    int i , RET;
    int n_variable_entries = 2;
    int* varying_cols = (int*)calloc(n_variable_entries, sizeof(int));
    int* varying_rows = (int*)calloc(n_variable_entries, sizeof(int));
    varying_cols[0] = 3;
    varying_cols[1] = 8;
    varying_rows[0] = 4;
    varying_rows[1] = 6;

    klu_defaults (&Common) ;

    /* preprocess matrix */
    Symbolic = klu_analyze_partial (n, Ap, Ai, varying_cols, varying_rows, n_variable_entries, AMD_ORDERING_RA, &Common) ;
    if(!Symbolic)
    {
        goto FAIL;
    }

    /* numerically factor the matrix */
    Numeric = klu_factor (Ap, Ai, Ax, Symbolic, &Common) ;
    if(!Numeric)
    {
        goto FAIL;
    }

    /* determine first varying column */
    RET = klu_determine_start(Symbolic, Numeric, &Common, Ap, Ai, varying_cols, varying_rows, n_variable_entries);
    if(RET != (1))
    {
        goto FAIL;
    }
    printf("------------------\n");
    printf("------------------\n");
    
    /* partially refactor based on refactorization restart */
    RET = klu_partial_refactorization_restart(Ap, Ai, Ax_new, Symbolic, Numeric, &Common);
    if(RET != (1))
    {
        goto FAIL;
    }

    /* solve rhs */
    klu_solve (Symbolic, Numeric, 10, 1, b, &Common) ;
    printf("Partial solution:\n");
    for (i = 0 ; i < n ; i++) printf ("xp [%d] = %g\n", i, b [i]) ;
    printf("------------------\n");

    /* compute full refactorization */
    RET = klu_refactor(Ap, Ai, Ax_new, Symbolic, Numeric, &Common);
    if(RET != (1))
    {
        goto FAIL;
    }

    /* solve rhs */
    RET = klu_solve (Symbolic, Numeric, 10, 1, c, &Common) ;
    if(RET != (1))
    {
        goto FAIL;
    }

    printf("------------------\n");
    printf("Refactor solution:\n");
    for (i = 0 ; i < n ; i++) printf ("xr [%d] = %g\n", i, c [i]) ;
    printf("------------------\n");

    printf("------------------\n");
    printf("Difference:\n");
    for(i = 0; i < n ; i++) printf("xp[%d] - xr[%d] = %g\n", i, i, fabs(b[i]-c[i]));
    printf("------------------\n");

    /* check for errors */
    max_error = fabs(b[0]-c[0]);
    for(i = 1 ; i < n ; i++)
    {
        error = fabs(b[i]-c[i]);
        if(error > max_error)
        {
            max_error = error;
        }
    }

    if(error > TOLERANCE)
    {
        return (1);
    }

    klu_free_symbolic (&Symbolic, &Common) ;
    klu_free_numeric (&Numeric, &Common) ;
    free(varying_cols);
    free(varying_rows);
    return (0) ;

FAIL:
    klu_free_symbolic (&Symbolic, &Common) ;
    klu_free_numeric (&Numeric, &Common) ;
    free(varying_cols);
    free(varying_rows);
    return (1) ;
}

