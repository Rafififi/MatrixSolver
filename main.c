#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        printf("Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    // <Handle the inputs here>
    const char *filename = argv[1];

    CSRMatrix aMatrix;
    ReadMMtoCSR(filename, &aMatrix);

    // Initializing all the vector b (in Ax=b)
    double *bMatrix = (double *)malloc(aMatrix.num_cols * sizeof(double));
    double *xMatrix = (double *)malloc(aMatrix.num_cols * sizeof(double));
    // Set all elements of b to 1
    for (int i = 0; i < aMatrix.num_cols; ++i)
    {
        bMatrix[i] = 1.0;
    }

    for (int i = 0; i < aMatrix.num_cols; ++i)
    {
        xMatrix[i] = 1.0;
    }

    
    double *r = (double *)malloc(aMatrix.num_cols * sizeof(double));
    compute_residual(aMatrix, bMatrix, xMatrix, r);
    spmv_csr(&aMatrix, xMatrix, r);
    double norm = compute_norm(r, aMatrix.num_cols);
    printf("Norm: %lf\n", norm);
    // <The rest of your code goes here>
    free(bMatrix);
    free(xMatrix);
    free(r);

    free_csr_matrix(&aMatrix);

    printf("\nDone!\n");
    return 0;
}