#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include <string.h>
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
    clock_t start, midpoint1, midpoint2, end;
    double cpu_time_used;
    start = clock();


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

    CSRMatrix *aMatrixTranspose = (CSRMatrix *)malloc(sizeof(CSRMatrix)); // allocate memory for the matrix that will store the transpose of the original matrix

    // check if memory allocation was successful
    if (aMatrixTranspose == NULL)
    {
        printf("Error: memory allocation failed\n");
        exit(1);
    }

    char matrixType = triangularCheck(aMatrix, aMatrixTranspose);


   
    midpoint1 = clock(); // start timing the solver
    solver(aMatrix, bMatrix, xMatrix, *aMatrixTranspose);
    midpoint2 = clock(); // end timing the solver
    puts("Solver finished");

    if (aMatrix.num_non_zeros < 51)
    {
        for (int i = 0; i < aMatrix.num_cols; ++i)
        {
            printf("%e\n", xMatrix[i]);
        }
    }
    
    double *residual = (double *)malloc(aMatrix.num_cols * sizeof(double));
    if (residual == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }

    computeResidual(aMatrix, bMatrix, xMatrix, residual, *aMatrixTranspose);
    

    
    double norm = computeNorm(residual, aMatrix.num_cols);
    printf("Norm: %e\n", norm);
    // <The rest of your code goes here>
    freeCSRMatrix(&aMatrix);
    if (matrixType == 'L')
    {
        freeCSRMatrix(aMatrixTranspose);
    }
    free(aMatrixTranspose);
    free(bMatrix);
    free(xMatrix);
    free(residual);
    bMatrix = NULL;
    xMatrix = NULL;
    residual = NULL;

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken to run the program: %f\n", cpu_time_used);
    printf("Time taken to solve the matrix: %f\n", ((double)(midpoint2 - midpoint1)) / CLOCKS_PER_SEC);
    printf("\nDone!\n");
    return 0;
}