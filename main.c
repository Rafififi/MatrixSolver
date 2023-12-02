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
    spmv_csr(&aMatrix, xMatrix, r);
    // <The rest of your code goes here>
    free(bMatrix);
    free(xMatrix);
    free(r);
    
    printf("\nDone!\n");
    return 0;
}