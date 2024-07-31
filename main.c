#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <time.h>
#include <string.h>
#include "functions.h"


int main(int argc, char* argv[]) {

  if (argc != 2) {
    printf("Usage: %s [martix-market-filename]\n", argv[0]);
    exit(1);
  }
  const char* filename = argv[1];
  const clock_t start = clock();

  CSRMatrix aMatrix;
  ReadMMtoCSR(filename, &aMatrix);

  puts("Matrix read successfully");
  // Initializing all the vector b (in Ax=b)
  double* bMatrix = malloc(aMatrix.num_cols * sizeof(double));
  double* xMatrix = malloc(aMatrix.num_cols * sizeof(double));
  if (!bMatrix || !xMatrix){ memFail(); }

  // Set all elements of b to 1
  for (int i = 0; i < aMatrix.num_cols; i++) {
    bMatrix[i] = 1.0;
  }

  for (int i = 0; i < aMatrix.num_cols; i++) {
    xMatrix[i] = 1.0;
  }

  CSRMatrix* aMatrixTranspose = malloc(sizeof(CSRMatrix)); // allocate memory for the matrix that will store the transpose of the original matrix
  if (!aMatrixTranspose){ memFail(); }

  char matrixType = triangularCheck(aMatrix, aMatrixTranspose);

  const clock_t midpoint1 = clock();
  solver(aMatrix, bMatrix, xMatrix, *aMatrixTranspose);
  const clock_t midpoint2 = clock(); 
  puts("Solver finished");

  if (aMatrix.num_non_zeros < 51) {
    for (int i = 0; i < aMatrix.num_cols; ++i) {
      printf("%e\n", xMatrix[i]);
    }
  }

  double* residual = malloc(aMatrix.num_cols * sizeof(double));
  if (!residual){  memFail(); }

  computeResidual(aMatrix, bMatrix, xMatrix, residual, *aMatrixTranspose);

  double norm = computeNorm(residual, aMatrix.num_cols);
  printf("Norm: %e\n", norm);
  freeCSRMatrix(&aMatrix);
  if (matrixType == 'L') {
    freeCSRMatrix(aMatrixTranspose);
  }
  free(aMatrixTranspose);
  free(bMatrix);
  free(xMatrix);
  free(residual);
  bMatrix = NULL;
  xMatrix = NULL;
  residual = NULL;

  const clock_t end = clock();
  double cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
  printf("Time taken to run the program: %f\n", cpu_time_used);
  printf("Time taken to solve the matrix: %f\n", ((double)(midpoint2 - midpoint1)) / CLOCKS_PER_SEC);
  printf("\nDone!\n");
  return 0;

}
