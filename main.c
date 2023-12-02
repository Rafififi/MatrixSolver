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
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("File not found\n");
        exit(1);
    }

    CSRMatrix aMatrix;
    ReadMMtoCSR(filename, &aMatrix);
    puts("Matrix read successfully");
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

    //solver(aMatrix, bMatrix, xMatrix);
    double *r = (double *)malloc(aMatrix.num_cols * sizeof(double));
    compute_residual(aMatrix, bMatrix, xMatrix, r);
    double norm = compute_norm(r, aMatrix.num_cols);
    printf(" \nNorm: %lf\n", norm);
    // <The rest of your code goes here>
    free_csr_matrix(&aMatrix);

    free(bMatrix);
    free(xMatrix);
    free(r);
    bMatrix = NULL;
    xMatrix = NULL;
    r = NULL;

    printf("\nDone!\n");
    return 0;
}