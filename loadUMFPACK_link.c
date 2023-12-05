#include <stdio.h> //perror
#include <dlfcn.h>
#include "umfpack.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/*
gcc -o result loadUMFPACK_link.c -I/Users/walter/Documents/SuiteSparse/UMFPACK/Include -I/Users/walter/Documents/SuiteSparse/SuiteSparse_config -I/Users/walter/Documents/SuiteSparse/AMD/Include  -L
/Users/walter/github/SuiteSparse/stable/lib/ -lumfpack -lm -lpthread
*/
/*********************************************************************************
 * failed_matrix.dat is a binary that stores our instances.
 *  
 * failed_matrix.dat contains the matrix information in the 
 * following order :
 * 
 * 1. int size    (the size of the square matrix. i.e the number of rows and columns)
 * 2. int nnz     (number of non-zero element in our data matrix)
 * 3. int *Ai     (An array of length 'nnz' to store the row indices
 * 4. int *Ap     (An array of length 'size + 1')
 * 5. double *Ax  (An array of length 'nnz' to store matrix data)
 * 
 * 
 * So it stores everything based on UMFPACK's form.
 * *******************************************************************************/


/*********************************************************************************
 * This file does the following steps 
 * 1. Read in the matrix that fails to be factorized.
 * 2. Dynamically load the routines for factorization. (you need to input your
 *    UMFPACK_PATH to link to the dynamic libraray)
 * 3. Try to factorize it using int32 interface.
 * 4. Failed with Out-Of-Memory Error
 ********************************************************************************/

int main(int argc, char **argv){


    FILE *fptr;
    int status;
    int size;
    int nnz;
    int *Ai;
    int *Ap;
    double *Ax;
    double *null = (double *) NULL ;
    void *Symbolic, *Numeric ;

    /* read in matrix */

    /* open file*/
    fptr = fopen("failed_matrix.dat","r");
    if (fptr == NULL){
        perror("file open failed");
    }

    /* read in all matrix information */
    fread(&size, sizeof(int), 1, fptr);
    fread(&nnz, sizeof(int), 1, fptr);

    Ai = (int *)malloc(nnz*sizeof(int));
    Ap = (int *)malloc((size + 1)*sizeof(int));
    Ax = (double *)malloc(nnz*sizeof(double));

    fread(Ai, sizeof(int), nnz, fptr);
    fread(Ap, sizeof(int), size + 1, fptr);
    fread(Ax, sizeof(double), nnz, fptr);
    /* finish reading matrix*/


    /* start factorization step*/

    printf("run symbolic step\n");
    status = umfpack_di_symbolic(size, size, Ap, Ai, Ax, &Symbolic, null, null);
    if (status != 0){
        perror("symbolic");
    }

    printf("run numeric step\n");
    status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
    if (status != 0){
        printf("status code : %d\n", status);
    }

    free(Ai);
    free(Ap);
    free(Ax);
    return 0;
}