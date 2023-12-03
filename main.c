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

    //If the file is not in the directory
    

    CSRMatrix aMatrix;
    ReadMMtoCSR(filename, &aMatrix);

    
    puts("Matrix read successfully");
    // Initializing all the vector b (in Ax=b)
    double *bMatrix = (double *)malloc(aMatrix.num_cols * sizeof(double));
    double *xMatrix = (double *)malloc(aMatrix.num_cols * sizeof(double));
    if (bMatrix == NULL || xMatrix == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }

    // Set all elements of b to 1
    for (int i = 0; i < aMatrix.num_cols; ++i)
    {
        bMatrix[i] = 1.0;
    }

    for (int i = 0; i < aMatrix.num_cols; ++i)
    {
        xMatrix[i] = 1.0;
    }


   

    //solver(aMatrix, bMatrix, xMatrix);
    puts("Solver finished");
    double *residual = (double *)malloc(aMatrix.num_cols * sizeof(double));
    if (residual == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }

    computeResidual(aMatrix, bMatrix, xMatrix, residual);
    double norm = computeNorm(residual, aMatrix.num_cols);
    printf(" \nNorm: %e\n", norm);
    // <The rest of your code goes here>
    freeCSRMatrix(&aMatrix);

    free(bMatrix);
    free(xMatrix);
    free(residual);
    bMatrix = NULL;
    xMatrix = NULL;
    residual = NULL;

    printf("\nDone!\n");
    return 0;
}