#include "klu_internal.h"
#include <string.h>
#define PRECISION 20

/* definitions for printing. if data type Int is defined to be long int, %ld printing is used. otherwise,
 * %d printing. */
#ifdef DLONG
    #define TYPE ("%ld")
    #define TYPEC ("%ld, ")
    #define TYPEN ("%ld\n")
    #define TYPEH ("%ld %ld %ld\n")
#ifdef COMPLEX
    #define TYPEF ("%ld %d %.*e %.*e\n")
    #define MATRIXMARKET_HEADER ("%%%%MatrixMarket matrix coordinate complex general\n")
#else
    #define TYPEF ("%ld %d %.*e\n")
    #define MATRIXMARKET_HEADER ("%%%%MatrixMarket matrix coordinate real general\n")
#endif
#else
    #define TYPE ("%d")
    #define TYPEC ("%d, ")
    #define TYPEN ("%d\n")
    #define TYPEH ("%d %d %d\n")
#ifdef COMPLEX
    #define TYPEF ("%d %d %.*e %.*e\n")
    #define MATRIXMARKET_HEADER ("%%%%MatrixMarket matrix coordinate complex general\n")
#else
    #define TYPEF ("%d %d %.*e\n")
    #define MATRIXMARKET_HEADER ("%%%%MatrixMarket matrix coordinate real general\n")
#endif
#endif

/* prints out BTF+AMD+PP (or any other ordering) permutation WITH partial pivoting */

void KLU_dumpPerm(Int* Q, Int* P, Int n, Int counter)
{
    int i;
    char strQ[32];
    char strP[32];
    char counterstring[32];
    sprintf(counterstring, TYPE, counter);

    strcpy(strQ, "KLU_Q");
    strcpy(strP, "KLU_P");
    strcat(strQ, counterstring);
    strcat(strP, counterstring);
    strcat(strQ, ".txt");
    strcat(strP, ".txt");

    /* dump permutations into file */
    FILE* fQ = fopen(strQ, "w");
    FILE* fP = fopen(strP, "w");

    /* dump Q */
    for (i = 0 ; i < n-1; i++)
    {
        fprintf(fQ, TYPEC, Q[i]);
    }
    fprintf(fQ, TYPEN, Q[n-1]);

    /* dump P */
    for (i = 0 ; i < n-1; i++)
    {
        fprintf(fP, TYPEC, P[i]);
    }
    fprintf(fP, TYPEN, P[n-1]);

    fclose(fQ);
    fclose(fP);
}

/* prints out BTF+AMD (or any other ordering) permutation WITHOUT partial pivoting */
void KLU_dumpPermPre(Int* Q, Int* P, Int n, Int counter)
{
    int i;
    char strQ[32];
    char strP[32];
    char counterstring[32];
    sprintf(counterstring, TYPE, counter);

    strcpy(strQ, "KLU_QPrev");
    strcpy(strP, "KLU_PPrev");
    strcat(strQ, counterstring);
    strcat(strP, counterstring);
    strcat(strQ, ".txt");
    strcat(strP, ".txt");

    /* dump permutations into file */
    FILE* fQ = fopen(strQ, "w");
    FILE* fP = fopen(strP, "w");

    /* dump Q */
    for (i = 0 ; i < n-1; i++)
    {
        fprintf(fQ, TYPEC, Q[i]);
    }
    fprintf(fQ, TYPEN, Q[n-1]);

    /* dump P */
    for (i = 0 ; i < n-1; i++)
    {
        fprintf(fP, TYPEC, P[i]);
    }
    fprintf(fP, TYPEN, P[n-1]);

    fclose(fQ);
    fclose(fP);
}

void KLU_printMTX(FILE* f, double* values, Int* indices, Int* pointers, Int n, Int nz)
{
    int i,j;
    fprintf(f, MATRIXMARKET_HEADER);
    fprintf(f, TYPEH, n, n, nz);
    for(i = 0; i < n ; i++)
    {
        for(j = pointers[i] ; j < pointers[i+1] ; j++)
        {
            fprintf(f, TYPEF, indices[j]+1, i+1, PRECISION, values[j]);
        }
    }
}

void KLU_dumpLU(double *Lx, Int *Li, Int *Lp, 
            double *Ux, Int *Ui, Int *Up, 
            double *Fx, Int *Fi, Int *Fp, 
            Int lnz, Int unz, Int n, 
            Int nzoff, Int counter)
{
    char strL[32];
    char strU[32];
    char strF[32];
    char counterstring[32];
    sprintf(counterstring, TYPE, counter);
    strcpy(strL, "KLU_L");
    strcpy(strU, "KLU_U");
    strcpy(strF, "KLU_F");
    strcat(strL, counterstring);
    strcat(strU, counterstring);
    strcat(strF, counterstring);
    strcat(strL, ".mtx");
    strcat(strU, ".mtx");
    strcat(strF, ".mtx");

    FILE *l, *u, *g;
    g = fopen(strF, "w");
    l = fopen(strL, "w");
    u = fopen(strU, "w");

    KLU_printMTX(g, Fx, Fi, Fp, n, nzoff);
    KLU_printMTX(l, Lx, Li, Lp, n, lnz);
    KLU_printMTX(u, Ux, Ui, Up, n, unz);

    fclose(g);
    fclose(l);
    fclose(u);
}

void KLU_dumpPath(Int* path, Int* bpath, Int n, Int nb, Int pathLen, Int counter)
{
        int i;
        char counterstring[32];
        sprintf(counterstring, TYPE, counter);
        char strbpath[32];
        char strpath[32];
        strcpy(strbpath, "KLU_bpath");
        strcpy(strpath, "KLU_path");
        strcat(strbpath, counterstring);
        strcat(strpath, counterstring);
        strcat(strbpath, ".txt");
        strcat(strpath, ".txt");
        FILE* fbpath = fopen(strbpath, "w");
        FILE* fpath = fopen(strpath, "w");

        /* dump path into file */
        for (i = 0; i < nb - 1; i++)
        {
            fprintf(fbpath, TYPEC, bpath[i]);
        }
        fprintf(fbpath, TYPEN, bpath[nb-1]);

        if(pathLen > 0)
        {
            for (i = 0; i < pathLen - 1; i++)
            {
                fprintf(fpath, TYPEC, path[i]);
            }
            fprintf(fpath, TYPEN, path[pathLen-1]);
        }

        fclose(fpath);
        fclose(fbpath);
}

void KLU_dumpAll(double *Lx, 
            Int *Li, 
            Int *Lp,
            double *Ux, 
            Int *Ui, 
            Int *Up, 
            double *Fx, 
            Int *Fi, 
            Int *Fp, 
            Int *P,
            Int *Q,
            Int *path,
            Int *bpath,
            Int lnz,
            Int unz,
            Int n,
            Int nzoff,
            Int nb,
            Int pathLen
        )
{
    static Int counter = 0;
    KLU_dumpPerm(Q, P, n, counter);
    KLU_dumpKLU(Lx, Li, Lp, Ux, Ui, Up, Fx, Fi, Fp, lnz, unz, n, nzoff, counter);
    KLU_dumpPath(path, bpath, n, nb, pathLen, counter);
    counter++;
}

void KLU_dumpA(double* Ax,
            Int* Ai,
            Int* Ap,
            Int n
        )
{
    static Int counter = 0;
    char str[32];
    char counterstring[32];
    sprintf(counterstring, TYPE, counter);
    strcpy(str, "KLU_A");
    strcat(str, counterstring);
    strcat(str, ".mtx");

    FILE *a;
    a = fopen(str, "w");

    KLU_printMTX(a, Ax, Ai, Ap, n, Ap[n]);

    fclose(a);
    counter++;
}
