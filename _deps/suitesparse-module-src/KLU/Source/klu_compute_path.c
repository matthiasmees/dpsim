/* ========================================================================== */
/* === KLU_compute_path ===================================================== */
/* ========================================================================== */

/*
 * Computes the factorization path given change vector
 * Any new refactorisation can be computed by iterating over entries in factorization path
 * instead of all columns (see klu_partial.c)
 */
#include "klu_internal.h"

Int KLU_compute_path(
                    KLU_symbolic *Symbolic,
                    KLU_numeric *Numeric,
                    KLU_common *Common,
                    Int Ap [ ],
                    Int Ai [ ],
                    Int *variable_columns,
                    Int *variable_rows,
                    Int n_variable_entries
                    )
{

    /*
     * Computes Factorization Path.
     * Factorization-Path mode partial refactorization works as follows:
     * five arrays for two cases are needed.
     * for a variable entry, compute the permutation, i.e. A to LU
     *                  - case 1: entry is in off-diagonal block F
     *                          - KLU handles these blocks annoyingly
     *                          - two new arrays are used for partial refactorization:
     *                              - variable_offdiag_orig_entry[i]: refers to the i-th entry in Ax that is varying
     *                              - variable_offdiag_perm_entry[i]: refers to the corresponding position of i in F
     *                          - these entries are excluded from factorization path computation
     *                          - in partial refactorization, these values are simply "copied" to the correct location
     *                              - this can then be done very efficiently, i.e. with parallelization or SIMD (maybe)
     *                  - case 2: entry is in a diagonal block
     *                          - three arrays are used:
     *                              - variable_block: has length n_variable_blocks, contains indices of variable blocks
     *                                  - e.g. variable_block = {0, 3} means that blocks 0 and 3 contain variable entries
     *                              - block_path: has length nblocks (number of blocks), indicates position of varying columns for a block in the factorization path
     *                                  - a bit more complicated
     *                                  - e.g. block_path = {0, 0, 3, 3} means that block 1's variable columns are from index 0 to 2 in factorization path (variable_block = {1})
     *                              - path: factorization path
     *                                  - self-explanatory
     *                                  - e.g. path = {3, 4}
     *                           - in total:
     *                                  - variable_block = {1}
     *                                  - block_path = {0, 0, 3, 3}
     *                                  - path = {3, 4, 5}
     *                                  - means that block 1 contains varying entries. These are from index 0 to 2 in factorization path, which evaluates to columns 3, 4 and 5
     */

    if(variable_columns == NULL)
    {
        return FALSE;
    }
    if(variable_rows == NULL)
    {
        return FALSE;
    }
    if(n_variable_entries <= 0)
    {
        return TRUE;
    }

    /* This function is very long, because cases with and without BTF are implemented */
    /* Declarations */
    /* LU data */

    int n = Symbolic->n;
    int lnz = Numeric->lnz;
    int unz = Numeric->unz;
    int nzoff = Numeric->nzoff;
    int nb = Symbolic->nblocks;
    Int *Pinv = Numeric->Pinv;
    Int *Lp, *Li, *Up, *Ui, *Fi, *Fp;
    double *Lx, *Ux, *Fx;
    Int *P, *Q, *R;
    double *Rs;
    int RET;
    int oldcol, pend, p, newrow;
    int poff = 0, ctr = 0, variable_offdiag_length = 0;

    Int n_variable_entries_new = n_variable_entries, n_variable_columns = 0;

    /* arrays for the offdiagonal entries. are empty if BTF is not used */
    Int* variable_offdiag_orig_entry = (Int*)calloc(nzoff, sizeof(Int));
    Int* variable_offdiag_perm_entry = (Int*)calloc(nzoff, sizeof(Int));
    int* workpath = (int*)calloc(n, sizeof(int));
    Int* path = (Int*)calloc(n, sizeof(Int)); 
    Int* nvblocks = (Int*)calloc(nb, sizeof(Int));

    /* indices and temporary variables */
    int i, k, j, block;
    int pivot;
    int col = 0, nextcol;
    int flag = 1;

    /* blocks */
    Int k2 = n, k1 = 0, nk = 1;

    Int *Qi = (Int*)calloc(n, sizeof(Int));
    if(Qi == NULL)
    {
        return KLU_OUT_OF_MEMORY;
    }

    /* variable columns and rows in A are given. We need to know those in LU */
    Int *variable_columns_in_LU = (Int*)calloc(n_variable_entries, sizeof(Int));
    if(variable_columns_in_LU == NULL)
    {
        return KLU_OUT_OF_MEMORY;
    }

    Int *variable_rows_in_LU = (Int*)calloc(n_variable_entries, sizeof(Int));
    if(variable_rows_in_LU == NULL)
    {
        return KLU_OUT_OF_MEMORY;
    }

    /* if computation of fact. path (i.e. this function) was done before -> free memory */
    if (Numeric->path)
    {
        KLU_free(Numeric->path, n, sizeof(int), Common);
    }
    if (Numeric->block_path)
    {
        KLU_free(Numeric->block_path, Numeric->nblocks, sizeof(int), Common);
    }
    if (Numeric->variable_block)
    {
        KLU_free(Numeric->variable_block, Numeric->n_variable_blocks, sizeof(int), Common);
    }
    if (Numeric->variable_offdiag_orig_entry)
    {
        KLU_free(Numeric->variable_offdiag_orig_entry, Numeric->variable_offdiag_length, sizeof(int), Common);
    }
    if (Numeric->variable_offdiag_perm_entry)
    {
        KLU_free(Numeric->variable_offdiag_perm_entry, Numeric->variable_offdiag_length, sizeof(int), Common);
    }

    /* block path can be allocated now since its size is known */
    Numeric->block_path = KLU_malloc(nb+1, sizeof(int), Common);

    for (i = 0; i < nb; i++)
    {
        Numeric->block_path[i] = 0;
    }

    Lp = (Int*) calloc(n + 1, sizeof(Int));
    Up = (Int*) calloc(n + 1, sizeof(Int));
    Fp = (Int*) calloc(n + 1, sizeof(Int));
    Lx = (double*) calloc(lnz, sizeof(double));
    Ux = (double*) calloc(unz, sizeof(double));
    Fx = (double*) calloc(nzoff, sizeof(double));
    Li = (Int*) calloc(lnz, sizeof(Int));
    Ui = (Int*) calloc(unz, sizeof(Int));
    Fi = (Int*) calloc(nzoff, sizeof(Int));
    P = (Int*) calloc(n, sizeof(Int));
    Q = (Int*) calloc(n, sizeof(Int));
    Rs = (double*) calloc(n, sizeof(double));
    R = (Int*) calloc(nb + 1, sizeof(Int));

    /* ---------------------------------------------------------------- */
    /* first, get LU decomposition */
    /* sloppy implementation, there might be a smarter way to do this */
    /* ---------------------------------------------------------------- */
    RET = KLU_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, Common);

    /* check if extraction of LU matrix broke */
    if (RET != (TRUE))
    {
        return (FALSE);
    }

    /* ---------------------------------------------------------------- */
    /* second, invert permutation vector */
    /* Q gives "oldcol", we need "newcol" */
    /* ---------------------------------------------------------------- */
    for (i = 0; i < n; i++)
    {
        Qi[Q[i]] = i;
    }

    /* ---------------------------------------------------------------- */
    /* third, apply permutation on variable_columns and rows */
    /* ---------------------------------------------------------------- */

    for (i = 0; i < n_variable_entries; i++)
    {
        variable_columns_in_LU[i] = Qi[variable_columns[i]];
        variable_rows_in_LU[i] = Pinv[variable_rows[i]];
    }

    /* ---------------------------------------------------------------- */
    /* fourth, determine variable off-diagonal entries */
    /* only necessary if BTF used */
    /* if BTF is not used, no off-diagonal blocks are available */
    /* ---------------------------------------------------------------- */

    /* maybe write this function externally, because it's needed again by pr-rr */
    if(Common->btf == TRUE)
    {
        /* iterate over all blocks */
        for(block = 0 ; block < nb ; block++)
        {
            k1 = R [block];
            k2 = R [block+1];
            nk = k2 - k1;
            for (k = 0 ; k < nk ; k++)
            {
                oldcol = Q [k+k1] ;
                pend = Ap [oldcol+1] ;
                for (p = Ap [oldcol] ; p < pend ; p++)
                {
                    newrow = Pinv [Ai [p]] - k1 ;
                    if (newrow < 0 && poff < nzoff)
                    {
                        /* entry in off-diagonal block */

                        /* check if entry is in variable column */
                        /* TODO */
                        for(i = 0 ; i < n_variable_entries ; i++)
                        {
                            if(variable_columns_in_LU[i] == k+k1 && variable_rows_in_LU[i] == Pinv [Ai [p]])
                            {
                                /* set to -1, because they're off-diagonal now */
                                /* indicates that they are "removed" from variable columns and row "list" */
                                /* subsequently, add to variable off-diag entries and reduce length of variable entries */
                                variable_columns_in_LU[i] = -1;
                                variable_rows_in_LU[i] = -1;

                                variable_offdiag_orig_entry[variable_offdiag_length] = p;
                                variable_offdiag_perm_entry[variable_offdiag_length++] = poff;

                                n_variable_entries_new--;
                                break;
                            }
                        }
                        poff++ ;
                    }
                }
            }
        }
    }

    /* safety first */
    ASSERT(variable_offdiag_length == n_variable_entries - n_variable_entries_new);

    Numeric->variable_offdiag_orig_entry = KLU_malloc(variable_offdiag_length, sizeof(Int), Common);
    Numeric->variable_offdiag_perm_entry = KLU_malloc(variable_offdiag_length, sizeof(Int), Common);
    Numeric->variable_offdiag_length = variable_offdiag_length;

    for(i = 0; i < variable_offdiag_length ; i++)
    {
        Numeric->variable_offdiag_orig_entry[i] = variable_offdiag_orig_entry[i];
        Numeric->variable_offdiag_perm_entry[i] = variable_offdiag_perm_entry[i];
    }

    if(n_variable_entries_new == 0)
    {
        /* very good case, no varying columns in LU blocks */
        /* free memory and return */
        Numeric->pathLen = 0;
        Numeric->n_variable_blocks = 0;
        Numeric->path = KLU_malloc(0, sizeof(int), Common);
        Numeric->variable_block = KLU_malloc(Numeric->n_variable_blocks, sizeof(int), Common);
        goto EXIT;
    }

    int* cV = (int*) calloc(n, sizeof(int));
    if(cV == NULL)
    {
        return KLU_OUT_OF_MEMORY;
    }

    for(i = 0; i < n_variable_entries ; i++)
    {
        if(variable_columns_in_LU[i] != -1)
        {
            cV[variable_columns_in_LU[i]] = 1;
        }
    }
    for( i = 0; i < n_variable_entries ; i++)
    {
        variable_columns_in_LU[i] = 0;
    }


    for ( i = 0 ; i < n ; i++)
    {
        if(cV[i] == 1)
        {
            n_variable_columns++;
        }
    }

    variable_columns_in_LU  = (Int*)realloc(variable_columns_in_LU, sizeof(Int)*n_variable_columns);
    if(variable_columns_in_LU == NULL)
    {
        return KLU_OUT_OF_MEMORY;
    }
    for ( i = 0; i< n_variable_columns ; i++)
    {
        variable_columns_in_LU[i] = 0;
    }
    for ( i = 0 ; i < n ; i++)
    {
        if(cV[i] != 0)
        {
            variable_columns_in_LU[ctr++] = i;
        }
    }
    free(cV);

    ASSERT(ctr == variable_entries_new);

    /* ---------------------------------------------------------------- */
    /* fifth, compute factorization path of blocks */
    /* ---------------------------------------------------------------- */

    if (Common->btf == FALSE)
    {
        Numeric->variable_block = KLU_malloc(1, sizeof(int), Common);
        Numeric->block_path[0] = 0;
        Numeric->variable_block[0] = 0;
        Numeric->n_variable_blocks = 1;
        /* no blocks */
        for (i = 0; i < n_variable_columns; i++)
        {
            flag = 1;
            pivot = variable_columns_in_LU[i];
            if (path[pivot] == 1)
            {
                /* already computed pivot? */
                continue;
            }
            path[pivot] = 1;

            /* GP-based version */
            while (pivot < n && flag)
            {
                flag = 0;
                nextcol = Up[pivot + 1];

                /* scan Ui for values in row "pivot" */
                for (j = nextcol; j < unz; j++)
                {
                    if (Ui[j] == pivot)
                    {
                        /* column j is affected
                            * find column in which j is
                        */
                        for (k = 0; k < n; k++)
                        {
                            if (j >= Up[k] && j < Up[k + 1])
                            {
                                col = k;
                                break;
                            }
                        }
                        /* col is always well-defined */
                        path[col] = 1;
                        workpath[col] = 1;
                    }
                }
                for (j = pivot + 1; j < n; j++)
                {
                    if (workpath[j] == 1)
                    {
                        workpath[j] = 0;
                        pivot = j;
                        flag = 1;
                        break;
                    }
                }
            }
        }
        ctr = 0;
        for( i = 0 ; i < n ; i++)
        {
            if(path[i] == 1)
            {
                ctr++;
            }
        }
        Numeric->pathLen = ctr;
        Numeric->path = KLU_malloc(ctr, sizeof(int), Common);
        for(i = 0; i < ctr ; i++)
        {
            Numeric->path[i] = 0;
        }
        Numeric->block_path[1] = ctr;
        ctr = 0;
        for(i = 0; i < n ; i++)
        {
            if(path[i] == 1)
            {
                Numeric->path[ctr++] = i;
            }
        }
    }
    else
    {
        Numeric->n_variable_blocks = 0;
        for (i = 0; i < n_variable_columns; i++)
        {
            flag = 1;
            /* get next changing column */
            pivot = variable_columns_in_LU[i];

            /* check if it was already computed */
            if (path[pivot] == 1)
            {
                /* already computed pivot
                 * do nothing, go to next pivot
                 */
                continue;
            }

            /* set first value of singleton path */
            path[pivot] = 1;

            /* find block of pivot */
            for (k = 0; k < nb; k++)
            {
                k1 = R[k];
                k2 = R[k + 1];
                if (k1 <= pivot && pivot < k2)
                {
                    nk = k2 - k1;

                    /* set varying block */

                    if(nvblocks[k] != 1)
                    {
                        nvblocks[k] = 1;
                        Numeric->n_variable_blocks += 1;
                    }
                    break;
                }
            }

            if (nk == 1)
            {
                /* 1x1-block, its pivot already in path */
                continue;
            }

            /* propagate until end
             * in blocks, pivot < n_block[k]
             *
             * k2 is the first column of the NEXT block
             *
             * if pivot == k2 - 1, it is the last column of the
             * current block. then, "nextcol" would be k2, which is already
             * in the next block. Thus, only do loop if pivot < k2 - 1
             */
            while (pivot < k2 - 1 && flag == 1)
            {
                flag = 0;
                nextcol = Up[pivot + 1];
                /* find closest off-diagonal entries in U */
                /* u saved in csc:
                   Ux = [x, x, x, x, x] entries saved column-wise
                   Ui = [0, 2, 1, 4, 3] row position of those entries in columns
                   Up = [0, 1, 2, 3, 4, 5] column pointers
                   to find closest off-diagonal of row k:
                   find all Ui[j] == pivot
                   find Up[k] <= j < Up[k+1], k > pivot
                   ALTERNATIVELY
                   Transpose U
                 */

                /* This loop finds all nnz in row "pivot",
                 * meaning their respective columns depend on pivot column
                 */
                for (j = nextcol; j < unz; j++)
                {
                    if (Ui[j] == pivot)
                    {
                        /* find column in which j is */
                        for (k = k1; k < k2; k++)
                        {
                            if (j >= Up[k] && j < Up[k + 1])
                            {
                                col = k;
                                path[col] = 1;
                                workpath[col] = 1;
                                break;
                            }
                        }
                    }
                }
                for (j = pivot + 1; j < k2; j++)
                {
                    if (workpath[j] == 1)
                    {
                        workpath[j] = 0;
                        pivot = j;
                        flag = 1;
                        break;
                    }
                }
            }
        }

        /* count how many variable columns there are */
        ctr = 0;
        for( i = 0 ; i < n ; i++)
        {
            if(path[i] == 1)
            {
                ctr++;
            }
        }
        Numeric->pathLen = ctr;
        /* allocate and initialize memory */
        Numeric->path = KLU_malloc(ctr, sizeof(int), Common);
        Numeric->variable_block = KLU_malloc(Numeric->n_variable_blocks, sizeof(int), Common);
        for(i = 0; i < ctr ; i++)
        {
            Numeric->path[i] = 0;
        }

        /* set variable blocks */
        ctr = 0;
        for(i = 0 ; i < nb ; i++)
        {
            if(nvblocks[i] == 1)
            {
                Numeric->variable_block[ctr++] = i;
            }
        }

        ASSERT(ctr == Numeric->n_variable_blocks);

        /* assemble factorization path */
        ctr = 0;
        for(i = 0; i < n ; i++)
        {
            if(path[i] == 1)
            {
                Numeric->path[ctr++] = i;
            }
        }
        block = 0;

        /* determine from where to where in each block there are variable columns */
        k = 0;
        for(i = 0 ; i < Numeric->n_variable_blocks ; i++)
        {
            k1 = R[Numeric->variable_block[i]];
            k2 = R[Numeric->variable_block[i]+1];
            Numeric->block_path[Numeric->variable_block[i]] = k;
            while(k < ctr && Numeric->path[k] < k2 - 1)
            {
                k++;
            }
            k++;
            Numeric->block_path[Numeric->variable_block[i]+1] = k;
        }
    }
EXIT:
    free(Lp);
    free(Li);
    free(Lx);
    free(Up);
    free(Ui);
    free(Ux);
    free(Fi);
    free(Fp);
    free(Fx);
    free(P);
    free(Q);
    free(Qi);
    free(R);
    free(Rs);
    free(variable_columns_in_LU);
    free(variable_rows_in_LU);
    free(variable_offdiag_orig_entry);
    free(variable_offdiag_perm_entry);
    free(workpath);
    free(path);
    free(nvblocks);
    return (TRUE);
}

/*
 * This function determines the first varying column for
 * partial refactorization by refactorization restart (PR-RR)
 */
Int KLU_determine_start(
        KLU_symbolic *Symbolic,
        KLU_numeric *Numeric,
        KLU_common *Common,
        Int Ap [ ],
        Int Ai [ ],
        Int *variable_columns,
        Int *variable_rows,
        Int n_variable_entries
    )
{
    int n = Symbolic->n;
    int lnz = Numeric->lnz;
    int unz = Numeric->unz;
    int nzoff = Numeric->nzoff;
    int nb = Symbolic->nblocks;
    Int *Pinv = Numeric->Pinv;
    Int *Lp, *Li, *Up, *Ui, *Fi, *Fp;
    double *Lx, *Ux, *Fx;
    Int *P, *Q, *R;
    double *Rs;
    int RET;

    int *variable_offdiag_orig_entry = (int*)calloc(nzoff, sizeof(int));
    int *variable_offdiag_perm_entry = (int*)calloc(nzoff, sizeof(int));
    int *nvblocks = (int*)calloc(nb, sizeof(int));

    Numeric->pathLen = 0;

    /* indices and temporary variables */
    int i, k, p, block;
    int pivot;
    int oldcol, pend, newrow, poff = 0, variable_offdiag_length = 0;
    int n_variable_entries_new = n_variable_entries;

    /* blocks */
    Int k2, k1, nk;

    Int *Qi = (Int*)calloc(n, sizeof(Int));
    if(Qi == NULL)
    {
        return KLU_OUT_OF_MEMORY;
    }

    Int *variable_columns_in_LU = (Int*)calloc(n_variable_entries, sizeof(Int));
    if(variable_columns_in_LU == NULL)
    {
        free(Qi);
        return KLU_OUT_OF_MEMORY;
    }

    Int *variable_rows_in_LU = (Int*)calloc(n_variable_entries, sizeof(Int));
    if(variable_rows_in_LU == NULL)
    {
        free(Qi);
        free(variable_columns_in_LU);
        return KLU_OUT_OF_MEMORY;
    }

    if (Numeric->path)
    {
        KLU_free(Numeric->path, n, sizeof(int), Common);
    }
    if (Numeric->block_path)
    {
        KLU_free(Numeric->block_path, Numeric->nblocks, sizeof(int), Common);
    }

    Numeric->block_path = KLU_malloc(nb+1, sizeof(int), Common);

    for (i = 0; i < nb; i++)
    {
        Numeric->block_path[i] = n;
    }

    Lp = (Int*) calloc(n + 1, sizeof(Int));
    Up = (Int*) calloc(n + 1, sizeof(Int));
    Fp = (Int*) calloc(n + 1, sizeof(Int));
    Lx = (double*) calloc(lnz, sizeof(double));
    Ux = (double*) calloc(unz, sizeof(double));
    Fx = (double*) calloc(nzoff, sizeof(double));
    Li = (Int*) calloc(lnz, sizeof(Int));
    Ui = (Int*) calloc(unz, sizeof(Int));
    Fi = (Int*) calloc(nzoff, sizeof(Int));
    P = (Int*) calloc(n, sizeof(Int));
    Q = (Int*) calloc(n, sizeof(Int));
    Rs = (double*) calloc(n, sizeof(double));
    R = (Int*) calloc(nb + 1, sizeof(Int));

    /* first, get LU decomposition
     * sloppy implementation, as there might be a smarter way to do this
     */
    RET = KLU_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, Common);

    /* check if extraction of LU matrix broke */
    if (RET != (TRUE))
    {
        free(Qi);
        free(variable_columns_in_LU);
        free(variable_rows_in_LU);
        return (FALSE);
    }

    /* second, invert permutation vector
     * Q gives "oldcol", we need "newcol"
     */
    for (i = 0; i < n; i++)
    {
        Qi[Q[i]] = i;
    }

    /* ---------------------------------------------------------------- */
    /* third, apply permutation on variable_columns and rows */
    /* ---------------------------------------------------------------- */

    for (i = 0; i < n_variable_entries; i++)
    {
        variable_columns_in_LU[i] = Qi[variable_columns[i]];
        variable_rows_in_LU[i] = Pinv[variable_rows[i]];
    }

    /* ---------------------------------------------------------------- */
    /* fourth, determine variable off-diagonal entries */
    /* only necessary if BTF used */
    /* if BTF is not used, no off-diagonal blocks are available */
    /* ---------------------------------------------------------------- */

    if(Common->btf == TRUE)
    {
        /* iterate over all blocks */
        for(block = 0 ; block < nb ; block++)
        {
            k1 = R [block];
            k2 = R [block+1];
            nk = k2 - k1;
            for (k = 0 ; k < nk ; k++)
            {
                oldcol = Q [k+k1] ;
                pend = Ap [oldcol+1] ;
                for (p = Ap [oldcol] ; p < pend ; p++)
                {
                    newrow = Pinv [Ai [p]] - k1 ;
                    if (newrow < 0 && poff < nzoff)
                    {
                        /* entry in off-diagonal block */

                        /* check if entry is in variable column */
                        /* TODO */
                        for(i = 0 ; i < n_variable_entries ; i++)
                        {
                            if(variable_columns_in_LU[i] == k+k1 && variable_rows_in_LU[i] == Pinv [Ai [p]])
                            {
                                /* set to -1, because they're off-diagonal now */
                                /* indicates that they are "removed" from variable columns and row "list" */
                                /* subsequently, add to variable off-diag entries and reduce length of variable entries */
                                variable_columns_in_LU[i] = -1;
                                variable_rows_in_LU[i] = -1;

                                variable_offdiag_orig_entry[variable_offdiag_length] = p;
                                variable_offdiag_perm_entry[variable_offdiag_length++] = poff;

                                n_variable_entries_new--;
                                break;
                            }
                        }
                        poff++ ;
                    }
                }
            }
        }
    }

    /* safety first */
    ASSERT(variable_offdiag_length == n_variable_entries - n_variable_entries_new);

    Numeric->variable_offdiag_orig_entry = KLU_malloc(variable_offdiag_length, sizeof(Int), Common);
    Numeric->variable_offdiag_perm_entry = KLU_malloc(variable_offdiag_length, sizeof(Int), Common);
    Numeric->variable_offdiag_length = variable_offdiag_length;

    for(i = 0; i < variable_offdiag_length ; i++)
    {
        Numeric->variable_offdiag_orig_entry[i] = variable_offdiag_orig_entry[i];
        Numeric->variable_offdiag_perm_entry[i] = variable_offdiag_perm_entry[i];
    }

    if(n_variable_entries_new == 0)
    {
        /* very good case, no varying columns in LU blocks */
        /* free memory and return */
        goto EXIT;
    }

    int* cV = (int*) calloc(n, sizeof(int));
    for ( i = 0 ; i < n_variable_entries ; i++)
    {
        if(variable_columns_in_LU[i] != -1)
        {
            cV[variable_columns_in_LU[i]] = 1;
        }
    }

    /* count how many variable columns there are (not necessarily equal to n_variable_entries_new) */
    n_variable_entries = 0;
    for ( i = 0 ; i < n ; i++)
    {
        if(cV[i] == 1)
        {
            n_variable_entries++;
        }
    }
    int ctr = 0;
    variable_columns_in_LU = (Int*) realloc(variable_columns_in_LU, sizeof(Int)*n_variable_entries);
    if(variable_columns_in_LU == NULL)
    {
        return KLU_OUT_OF_MEMORY;
    }
    for ( i = 0; i< n_variable_entries ; i++)
    {
        variable_columns_in_LU[i] = 0;
    }
    for ( i = 0 ; i < n ; i++)
    {
        if(cV[i] == 1)
        {
            variable_columns_in_LU[ctr] = i;
            ctr++;
        }
    }
    free(cV);

    if (Common->btf == FALSE)
    {
        /* no btf case
         * => identify first varying entry in entire matrix */
        Numeric->n_variable_blocks = 1;
        Numeric->variable_block = KLU_malloc(Numeric->n_variable_blocks, sizeof(int), Common);
        Numeric->variable_block[0] = 0;

        /* find minimum in variable_columns_in_LU */
        pivot = variable_columns_in_LU[0];
        for(i = 1; i < n_variable_entries ; i++)
        {
            if(pivot > variable_columns_in_LU[i])
            {
                pivot = variable_columns_in_LU[i];
            }
        }

        /* set first varying column internally */
        Numeric->block_path[0] = pivot;
    }
    else
    {
        Numeric->n_variable_blocks = 0;
        /* btf case
         *
         * iterate over variable_columns_in_LU
         * for each variable column
         * find its block
         *      put block in block-path
         *      if block larger than 1
         *          start[block] = min(start[block], column), if start[block] != 0
         */

         for(i = 0; i < n_variable_entries ; i++)
         {
            /* grab variable column number i */
            pivot = variable_columns_in_LU[i];

            /* find block */
            for (k = 0; k < nb; k++)
            {
                k1 = R[k];
                k2 = R[k + 1];
                if (k1 <= pivot && pivot < k2)
                {
                    nk = k2 - k1;

                    /* set varying block */

                    if(nvblocks[k] != 1)
                    {
                        nvblocks[k] = 1;
                        Numeric->n_variable_blocks += 1;
                    }
                    break;
                }
            }

            if (nk == 1)
            {
                /* 1x1-block, nothing left to do */
                continue;
            }

            /* nk x nk block. check if pivot column is smaller (before) current minimum */
            if(Numeric->block_path[k] > pivot)
            {
                Numeric->block_path[k] = pivot;
            }
         }
        ctr = 0;
        Numeric->variable_block = KLU_malloc(Numeric->n_variable_blocks, sizeof(int), Common);
        for(i = 0 ; i < nb ; i++)
        {
            if(nvblocks[i] == 1)
            {
                /* block i is varying */
                Numeric->variable_block[ctr++] = i;
            }
        }
        ASSERT(ctr == Numeric->n_variable_blocks);
    }
EXIT:
    free(Lp);
    free(Li);
    free(Lx);
    free(Up);
    free(Ui);
    free(Ux);
    free(Fi);
    free(Fp);
    free(Fx);
    free(P);
    free(Q);
    free(Qi);
    free(R);
    free(Rs);
    free(variable_columns_in_LU);
    free(variable_offdiag_orig_entry);
    free(variable_offdiag_perm_entry);
    free(nvblocks);

    return (TRUE);
}
