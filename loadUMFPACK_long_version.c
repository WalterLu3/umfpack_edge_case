#include <stdio.h> //perror
#include "umfpack.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>


#define Int64 long long

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
 * 2. Convert the matrix into int64 form. 
 * 3. Try to factorize it using int32 interface.
 * 4. numeric step never ends
 ********************************************************************************/

int main(int argc, char **argv){


    FILE *fptr;
    int status;
    int size;
    int nnz;
    int *Ai;
    int *Ap;
    int element;
    double *Ax;
    double *null = (double *) NULL ;
    void *Symbolic, *Numeric ;

    Int64 size_64;
    Int64 *Ai_64;
    Int64 *Ap_64;

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

    /* create 64-bits container */
    size_64 = (Int64) size;
    Ai_64 = (Int64 *)malloc(nnz*sizeof(Int64));
    Ap_64 = (Int64 *)malloc((size + 1)*sizeof(Int64));


    /* copy into 64-bits form */
    for (element = 1; element <= nnz; element++){
        Ai_64[element - 1] = (Int64) Ai[element - 1];
    }

    for (element = 1; element <= size + 1; element++){
        Ap_64[element - 1] = (Int64) Ap[element - 1];
    }
    

    /* start factorization step*/

    printf("run symbolic step\n");
    status = umfpack_dl_symbolic(size_64, size_64, Ap_64, Ai_64, Ax, &Symbolic, null, null);
    if (status != 0){
        perror("symbolic");
    }

    printf("run numeric step\n");
    status = umfpack_dl_numeric(Ap_64, Ai_64, Ax, Symbolic, &Numeric, null, null);
    if (status != 0){
        printf("status code : %d\n", status);
    }

    free(Ai);
    free(Ap);
    free(Ax);
    free(Ai_64);
    free(Ap_64);
    return 0;
}