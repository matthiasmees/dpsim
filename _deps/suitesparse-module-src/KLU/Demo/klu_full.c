/* klu_simpler: a simple KLU demo; solution is x = (1,2,3,4,5)xxx */

#include <stdio.h>
#include "klu.h"

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

void showLU(klu_symbolic* symb, klu_numeric* num, klu_common* com){
    int n = symb->n;
    int lnz = num->lnz;
    int unz = num->unz;
    int nzoff = num->nzoff;
    int *Lp, *Li, *Up, *Ui, *Fi, *Fp;
    double *Lx, *Ux, *Fx;
    int *P, *Q, *R;
    double* Rs;
    int nb = symb->nblocks;
    int i;
    Lp = calloc(n+1, sizeof(int));
    Up = calloc(n+1, sizeof(int));
    Fp = calloc(n+1, sizeof(int));
    Lx = calloc(lnz, sizeof(double));
    Ux = calloc(unz, sizeof(double));
    Fx = calloc(nzoff, sizeof(double));
    Li = calloc(lnz, sizeof(int));
    Ui = calloc(unz, sizeof(int));
    Fi = calloc(nzoff, sizeof(int));
    P = calloc(n, sizeof(int));
    Q = calloc(n, sizeof(int));
    Rs = calloc(n, sizeof(double));
    R = calloc(nb+1, sizeof(int));

    int RET = klu_extract(num, symb, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, com);
    printf("KLU extract: %d\n", RET);

    printf("--------------\nU:");
    printf(" unz: %d\n", unz);
    for(i=0; i<unz; i++)
        printf("%lf,\t", Ux[i]);
    printf("\n");
    for(i=0; i<unz; i++)
        printf("%d,\t", Ui[i]);
    printf("\n");
    for(i=0; i<n+1; i++)
        printf("%d,\t", Up[i]);
    printf("\n");

    printf("--------------\nL:");
    printf(" lnz: %d\n", lnz);
    for(i=0; i<lnz; i++)
        printf("%lf,\t", Lx[i]);
    printf("\n");
    for(i=0; i<lnz; i++)
        printf("%d,\t", Li[i]);
    printf("\n");
    for(i=0; i<n+1; i++)
        printf("%d,\t", Lp[i]);
    printf("\n");

    printf("--------------\nF:");
    printf(" nzoff: %d\n", nzoff);
    for(i=0; i<nzoff; i++)
        printf("%lf,\t", Fx[i]);
    printf("\n");
    for(i=0; i<nzoff; i++)
        printf("%d,\t", Fi[i]);
    printf("\n");
    for(i=0; i<n+1; i++)
        printf("%d,\t", Fp[i]);
    printf("\n");

    
    printf("--------------\nRs:\n");
    for(i=0; i<n; i++)
        printf("%lf,\t", Rs[i]);
    printf("\n");
    printf("--------------\nP:\n");
    for(i=0; i<n; i++)
        printf("%d,\t", P[i]);
    printf("\n");
    printf("--------------\nQ:\n");
    for(i=0; i<n; i++)
        printf("%d,\t", Q[i]);
    printf("\n");

    printf("--------------\nR:");
    printf(" nblocks: %d\n", symb->nblocks);
    for(i=0; i<symb->nblocks+1; i++)
        printf("%d\t", R[i]);
    printf("\n");

    free(Lp);
    free(Up);
    free(Fp);
    free(Li);
    free(Ui);
    free(Fi);
    free(Lx);
    free(Ux);
    free(Fx);
    free(P);
    free(Q);
    free(Rs);
    free(R);
}

int main (void)
{
    klu_symbolic *Symbolic ;
    klu_numeric *Numeric ;
    klu_common Common ;
    int i , RET;
    klu_defaults (&Common) ;

    Symbolic = klu_analyze (n, Ap, Ai, &Common) ;
    Numeric = klu_factor (Ap, Ai, Ax, Symbolic, &Common) ;
    
    RET = klu_full_partial(Ap, Ai, Ax_new, Symbolic, 3, Numeric, &Common);

    klu_solve (Symbolic, Numeric, 10, 1, b, &Common) ;
    printf("------------------\n");
    printf("Partial solution:\n");
    for (i = 0 ; i < n ; i++) printf ("xp [%d] = %g\n", i, b [i]) ;
    printf("------------------\n");

    RET = klu_refactor(Ap, Ai, Ax_new, Symbolic, Numeric, &Common);
    klu_solve (Symbolic, Numeric, 10, 1, c, &Common) ;
    printf("------------------\n");
    printf("Refactor solution:\n");
    for (i = 0 ; i < n ; i++) printf ("xr [%d] = %g\n", i, c [i]) ;
    printf("------------------\n");

    printf("------------------\n");
    printf("Difference:\n");
    for(i = 0; i < n ; i++) printf("xp[%d] - xr[%d] = %g\n", i, i, b[i]-c[i]);
    printf("------------------\n");
    klu_free_symbolic (&Symbolic, &Common) ;
    klu_free_numeric (&Numeric, &Common) ;
    return (0) ;
}

