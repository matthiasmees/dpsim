/* klu_simpler: a simple KLU demo; solution is x = (1,2,3,4,5)xxx */

#include <stdio.h>
#include "klu.h"

int    n = 3 ;
int    Ap [ ] = {0, 3, 6, 9} ;
int    Ai [ ] = { 0,  1,  2,  0,  1,  2,  0, 1, 2} ;
double Ax [ ] = {1.250, -1.0, -0.25, -1.00, 1.83333, -0.33333, -0.25, -0.33333, 1.58333} ;
double Ax_new [ ] = {1.250, -1.0, -0.25, -1.00, 1.83333, -0.33333, -0.25, -0.33333, 1.08333} ;
double b [ ] = {1.0, 0.0, 0.0} ;

int main (void)
{
    klu_symbolic *Symbolic ;
    klu_numeric *Numeric ;
    klu_common Common ;
    int i ;
    klu_defaults (&Common) ;
    //Common.scale = -1;
    Symbolic = klu_analyze (n, Ap, Ai, &Common) ;
    Numeric = klu_factor (Ap, Ai, Ax, Symbolic, &Common) ;
    //int ret = klu_partial(Ap, Ai, Ax, Symbolic, Numeric, &Common) ;
    int n = Symbolic->n;
    int lnz = Numeric->lnz;
    int unz = Numeric->unz;
    int nzoff = Numeric->nzoff;
    int *Lp, *Li, *Up, *Ui, *Fi, *Fp;
    double *Lx, *Ux, *Fx;
    int *P, *Q, *R;
    double* Rs;
    int nb = Symbolic->nblocks;

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

    int RET = klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, &Common);
    printf("KLU extract: %d\n", RET);

    for(int i=0; i<unz; i++)
        printf("%lf,\t", Ux[i]);
    printf("\n");
    for(int i=0; i<unz; i++)
        printf("%d,\t", Ui[i]);
    printf("\n");
    for(int i=0; i<n+1; i++)
        printf("%d,\t", Up[i]);
    printf("\n");

    for(int i=0; i<lnz; i++)
        printf("%lf,\t", Lx[i]);
    printf("\n");
    for(int i=0; i<lnz; i++)
        printf("%d,\t", Li[i]);
    printf("\n");
    for(int i=0; i<n+1; i++)
        printf("%d,\t", Lp[i]);
    printf("\n");
    
    for(int i=0; i<n; i++)
        printf("%lf,\t", Rs[i]);
    printf("\n");
    for(int i=0; i<n; i++)
        printf("%d,\t", P[i]);
    printf("\n");
    for(int i=0; i<n; i++)
        printf("%d,\t", Q[i]);
    printf("\n");

    for(int i=0; i<Symbolic->nblocks+1; i++)
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

    int* changeVector = (int*) calloc(sizeof(int), 1);
    int changeLen = 1;
    changeVector[0] = 0;
    klu_compute_path(Symbolic, Numeric, &Common, changeVector, changeLen);
    free(changeVector);
    // last column changed:
    //int PATH[] = {2};
    //int pathLen = 1;
    //RET = klu_partial(Ap, Ai, Ax_new, Symbolic, /*PATH, pathLen,*/ Numeric, &Common);
    //printf("KLU partial: %d\n", RET);
    //RET = klu_refactor(Ap, Ai, Ax_new, Symbolic, Numeric, &Common);
    //printf("KLU refact: %d\n", RET);

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

    RET = klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, &Common);
    printf("KLU extract: %d\n", RET);

    for(int i=0; i<unz; i++)
        printf("%lf,\t", Ux[i]);
    printf("\n");
    for(int i=0; i<unz; i++)
        printf("%d,\t", Ui[i]);
    printf("\n");
    for(int i=0; i<n+1; i++)
        printf("%d,\t", Up[i]);
    printf("\n");

    for(int i=0; i<lnz; i++)
        printf("%lf,\t", Lx[i]);
    printf("\n");
    for(int i=0; i<lnz; i++)
        printf("%d,\t", Li[i]);
    printf("\n");
    for(int i=0; i<n+1; i++)
        printf("%d,\t", Lp[i]);
    printf("\n");
    
    for(int i=0; i<n; i++)
        printf("%lf,\t", Rs[i]);
    printf("\n");
    for(int i=0; i<n; i++)
        printf("%d,\t", P[i]);
    printf("\n");
    for(int i=0; i<n; i++)
        printf("%d,\t", Q[i]);
    printf("\n");

    for(int i=0; i<Symbolic->nblocks+1; i++)
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

    klu_free_symbolic (&Symbolic, &Common) ;
    klu_free_numeric (&Numeric, &Common) ;
    //for (i = 0 ; i < n ; i++) printf ("x [%d] = %g\n", i, b [i]) ;
    return (0) ;
}

