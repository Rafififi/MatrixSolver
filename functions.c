#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "functions.h"

char matrixType = ' ';

void ReadMMtoCSR(const char *filename, CSRMatrix *aMatrix)
{
    FILE *fileName = fopen(filename, "r");
    if (fileName == NULL)
    {
        printf("File not found\n");
        exit(1);
    }

    fscanf(fileName, "%*[^\n]\n"); // Skip the first line
    //skip lines that start with a %
    char c = getc(fileName);
    while (c == '%')
    {
        fscanf(fileName, "%*[^\n]\n"); // Skip the line
        c = getc(fileName); // Get the next character
    }
    ungetc(c, fileName); // Goes back 1 character because the while loop goes 1 character too far

    // Read the number of rows, columns and non-zero elements
    fscanf(fileName, "%d %d %d", &aMatrix->num_rows, &aMatrix->num_cols, &aMatrix->num_non_zeros);

    printf("Number of rows: %d\n", aMatrix->num_rows); //For some reason vtune thinks that this line is allocating memory and not deallocating it
    //I checked with valgrind and all memory is being freed
    printf("Number of columns: %d\n", aMatrix->num_cols);
    printf("Number of non-zero elements: %d\n", aMatrix->num_non_zeros);
    
    
    //These are temporary pointers to store the rows, cols and the data accociated with it
    int *rows = (int *)calloc(aMatrix->num_non_zeros, sizeof(int));
    int *cols = (int *)calloc(aMatrix->num_non_zeros,  sizeof(int));
    double *data = (double *)calloc(aMatrix->num_non_zeros,  sizeof(double));


    for (int i = 0; i<aMatrix->num_non_zeros; i++)
    {
        fscanf(fileName, "%d %d %lf", &rows[i], &cols[i], &data[i]);
        c = getc(fileName);
    }

    printf("Numbers have been read\n");


    TemporaryElement *temporaryElement = (TemporaryElement *)calloc(aMatrix->num_non_zeros, sizeof(TemporaryElement));
    if (temporaryElement == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }
    /*
    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        printf("%d %d %lf\n", rows[i], cols[i], data[i]);
    }
    */

    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        temporaryElement[i].row = rows[i];
        temporaryElement[i].col = cols[i];
        temporaryElement[i].data = data[i];
    }

    //using qsort to sort the rows and applying the same permutation to cols and data
    qsort(temporaryElement, aMatrix->num_non_zeros, sizeof(TemporaryElement), compare);
    //qsort is a built in function that sorts the array elements in ascending order
    printf("qsort done\n");


    aMatrix->row_ptr = (int *)calloc((aMatrix->num_rows + 1), sizeof(int));
    aMatrix->col_ind = (int *)calloc(aMatrix->num_non_zeros, sizeof(int));
    aMatrix->csr_data = (double *)calloc(aMatrix->num_non_zeros, sizeof(double));

    if (aMatrix->row_ptr == NULL || aMatrix->col_ind == NULL || aMatrix->csr_data == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }

    int current_row = 0;
    int value = 0;
    int col = 0;

    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        while (current_row < temporaryElement[i].row && current_row <= aMatrix->num_rows)
        {
            aMatrix->row_ptr[current_row] = value;
            current_row++;
        }
        value++; 
        aMatrix->col_ind[col] = temporaryElement[i].col - 1; //The mtx file indexes start at 1 so we need to subtract 1 to get the correct index
        aMatrix->csr_data[col] = temporaryElement[i].data;
        col++;
    }

    // Fill the rest of row_ptr with the total number of non-zero elements
    while (current_row <= aMatrix->num_rows)
    {
        aMatrix->row_ptr[current_row] = value;
        current_row++;
    }
    
    
    free(rows);
    free(cols);
    free(data);
    free(temporaryElement);
    temporaryElement = NULL;
    rows = NULL;
    cols = NULL;
    data = NULL;
    fclose(fileName);

}

void solver(const CSRMatrix AMatrix, double *b, double *x, const CSRMatrix nonConsantMatrix)
{
    double *diagonal = (double *)calloc(AMatrix.num_rows, sizeof(double));

    //The larger the file the larger the amount of itterations needed to get a good approximation
    int maxIterations = 1000;

    printf("Generations: %d\n", maxIterations);

    if (diagonal == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < AMatrix.num_rows; i++) // Find the diagonal elements
    {
        for (int j = AMatrix.row_ptr[i]; j < AMatrix.row_ptr[i + 1]; j++)
        {
            if (AMatrix.col_ind[j] == i)
            {
                diagonal[i] = AMatrix.csr_data[j];
                if (diagonal[i] == 0)
                {
                    printf("The matrix is singular\n");
                    free(diagonal);
                    return;
                }
            }
        }
    }

    if (matrixType == 'L')
    {
        for (int iteration = 0; iteration < maxIterations; iteration++) // Go through the iterations
        {
            for (int j = 0; j < AMatrix.num_rows; j++) // Go through the rows
            {
                double sum = 0.0;
                for (int k = AMatrix.row_ptr[j]; k < AMatrix.row_ptr[j + 1]; k++) // go up to value of the next row pointer
                {
                    if (AMatrix.col_ind[k] != j)
                    {
                        sum += AMatrix.csr_data[k] * x[AMatrix.col_ind[k]]; // find the A*x for the non diagonal elements and the 
                    }
                }
                
                for (int k = nonConsantMatrix.row_ptr[j]; k < nonConsantMatrix.row_ptr[j + 1]; k++)
                {
                    if (nonConsantMatrix.col_ind[k] != j)
                    {
                        sum += nonConsantMatrix.csr_data[k] * x[nonConsantMatrix.col_ind[k]]; //find the AT*x for the non diagonal elements
                    }
                }

                x[j] = (b[j] - sum) / diagonal[j];
            }
        }
        free(diagonal);
        diagonal = NULL;
        return;
        
    }

    

    for (int i = 0; i < maxIterations; i++)
    {
        for (int j = 0; j < AMatrix.num_rows; j++)
        {
            double sum = 0;
            for (int k = AMatrix.row_ptr[j]; k < AMatrix.row_ptr[j + 1]; k++)
            {
                if (AMatrix.col_ind[k] != j)
                {
                    sum += AMatrix.csr_data[k] * x[AMatrix.col_ind[k]];
                }
            }
            x[j] = (b[j] - sum) / diagonal[j];
        }
    }

    free(diagonal);
    diagonal = NULL;
}


void spmvCSR(const CSRMatrix *AMatrix, const double *x, double *y, const CSRMatrix *nonConsantMatrix)
{
    if (matrixType == 'L')
    {
        for (int row = 0; row < AMatrix->num_rows; row++)
        {
            y[row] = 0;
        }
        double *product = (double *)calloc(AMatrix->num_non_zeros, sizeof(double));
        for (int row = 0; row < AMatrix->num_rows; row++)
        {

            for (int j = AMatrix->row_ptr[row]; j < AMatrix->row_ptr[row + 1]; j++)
            {
                if (AMatrix->col_ind[j] != row)
                {
                    product[row] = AMatrix->csr_data[j] * x[AMatrix->col_ind[j]];
                    y[row] += product[row];
                }
            }
            
            // printf(" %lf", y[row]);
        }



        printf("\n");
        for (int row = 0; row < nonConsantMatrix->num_rows; row++)
        {
            for (int j = nonConsantMatrix->row_ptr[row]; j < nonConsantMatrix->row_ptr[row + 1]; j++)
            {
                if (nonConsantMatrix->col_ind[j] != row)
                {
                    product[row] = nonConsantMatrix->csr_data[j] * x[nonConsantMatrix->col_ind[j]];
                    y[row] += product[row];
                }
                else
                {
                    product[row] = 0;
                }
            }
        }

        for (int row = 0; row < nonConsantMatrix->num_rows; row++) 
        {
            for (int j = nonConsantMatrix->row_ptr[row]; j < nonConsantMatrix->row_ptr[row + 1]; j++)
            {
                if (nonConsantMatrix->col_ind[j] == row)
                {
                    product[row] = nonConsantMatrix->csr_data[j] * x[nonConsantMatrix->col_ind[j]];
                    y[row] += product[row];
                }
                else 
                {
                    product[row] = 0;
                }
            }
        }

        free(product);
    }

    else
    {
        double *product = (double *)calloc(AMatrix->num_non_zeros, sizeof(double));
        if (product == NULL)
        {
            printf("Memory allocation failed\n");
            exit(1);
        }
        //printf("Product:");
        for (int row = 0; row < AMatrix->num_rows; row++)
        {

            for (int j = AMatrix->row_ptr[row]; j < AMatrix->row_ptr[row + 1]; j++)
            {
                product[row] += AMatrix->csr_data[j] * x[AMatrix->col_ind[j]];
            }
            y[row] = product[row];
            printf(" %lf", y[row]);
        }
        free(product);
    } 
    printf("\nspmv done\n");
}

void computeResidual(const CSRMatrix AMatrix, const double *b, const double *x, double *r, const CSRMatrix nonConsantMatrix)
{
    double *product = (double *)calloc(AMatrix.num_non_zeros, sizeof(double));

    if (product == NULL)
    {
        printf("Memory allocation failed\n");
        exit(1);
    }
    spmvCSR(&AMatrix, x, product, &nonConsantMatrix);

    for (int i = 0; i < AMatrix.num_rows; i++)
    {
        r[i] = product[i] - b[i];
    }
    
    free(product);
    printf("residual done\n");

}

double computeNorm(const double *r, int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
    {
        norm += r[i] * r[i];
    }
    norm = sqrt(norm);
    printf("norm done\n");
    return norm;

}

void freeCSRMatrix(CSRMatrix *aMatrix)
{
    free(aMatrix->row_ptr);
    free(aMatrix->col_ind);
    free(aMatrix->csr_data);
    aMatrix->row_ptr = NULL; //set the pointers to NULL to avoid dangling pointers
    aMatrix->col_ind = NULL;
    aMatrix->csr_data = NULL;
    aMatrix->num_rows = 0; //set the values to 0 to avoid any errors
    aMatrix->num_cols = 0;
    aMatrix->num_non_zeros = 0;
    printf("free done\n");
}

int compare(const void *a, const void *b) //this is for qsort to sort the rows and applying the same permutation to cols and data
{
    TemporaryElement *elementA = (TemporaryElement *)a;
    TemporaryElement *elementB = (TemporaryElement *)b;
    return elementA->row - elementB->row;
    //It works by returning a negative, zero, or positive integer, 
    //depending on whether the first argument is less than, equal to, or greater than the second argument.
}

char triangularCheck(const CSRMatrix aMatrix, CSRMatrix *aMatrixTranspose) //this is needed to check what type of matrix is being solved, as if it lower triangular we need to transpose to make it symetrical
{
    int upper = 1;
    //This Checks if there is any values in the section above the diagonal
    for (int i = 0; i < aMatrix.num_rows; i++)
    {
        for (int j = aMatrix.row_ptr[i]; j < aMatrix.row_ptr[i + 1]; j++)
        {
            if (aMatrix.col_ind[j] < i)
            {
                upper = 0;
                break;
            }
        }
        if (upper == 0) // if there is a value above the diagonal then it breaks the loop ensuring that we don't go through the whole matrix
        {
            break;
        }
    }
    
    int lower = 1;
    // This Checks if there is any values in the section under the diagonal
    for (int i = 0; i < aMatrix.num_rows; i++)
    {
        for (int j = aMatrix.row_ptr[i]; j < aMatrix.row_ptr[i + 1]; j++)
        {
            if (aMatrix.col_ind[j] > i)
            {
                lower = 0;
                break;
            }
        }
        if (lower == 0) // if there is a value under the diagonal then it breaks the loop ensuring that we don't go through the whole matrix
        {
            break;
        }
    }

    if (upper == 1 && lower == 1)
    {
        printf("The matrix not triangular\n");
        matrixType = 'N';
        return 'N';
    }
    else if (upper == 1)
    {
        printf("The matrix is upper triangular\n");
        matrixType = 'U';
        return 'U';
        
    }
    else if (lower == 1)
    {
        printf("The matrix is lower triangular\n");
        matrixType = 'L';
        CreateNonConstantMatrix(aMatrix, aMatrixTranspose);
        return 'L';
    }
    return 'N';
}

void CSRTranspose(CSRMatrix *aMatrix) //this is used to transpose the matrix if it is lower triangular
{
    // create temporary array of elements
    TemporaryElement *elements = (TemporaryElement *)calloc(aMatrix->num_non_zeros, sizeof(TemporaryElement));

    // check if memory allocation was successful
    if (elements == NULL)
    {
        printf("Error: memory allocation failed\n");
        return;
    }

    // iterate through the matrix and store the elements in the temporary array
    // the row and column of each element are swapped
    int index = 0;
    for (int i = 0; i < aMatrix->num_rows; i++)
    {
        for (int j = aMatrix->row_ptr[i]; j < aMatrix->row_ptr[i + 1]; j++)
        {
            elements[index].row = aMatrix->col_ind[j];
            elements[index].col = i;
            elements[index].data = aMatrix->csr_data[j];
            index++;
        }
    }
    

    // sort the temporary array by row
    qsort(elements, aMatrix->num_non_zeros, sizeof(TemporaryElement), compare);
    int current_row = 0;
    int value = 0;
    int col = 0;

    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        while (current_row < elements[i].row && current_row <= aMatrix->num_rows)
        {
            aMatrix->row_ptr[current_row] = value;
            current_row++;
        }
        value++;
        aMatrix->col_ind[col] = elements[i].col;
        aMatrix->csr_data[col] = elements[i].data;
        col++;
    }

    // Fill the rest of row_ptr with the total number of non-zero elements
    while (current_row <= aMatrix->num_rows)
    {
        aMatrix->row_ptr[current_row] = value;
        current_row++;
    }

    //shift all values of row_ptr to the right by 1 remove the final value and set the first value to 0
    for (int i = aMatrix->num_rows; i > 0; i--)
    {
        aMatrix->row_ptr[i] = aMatrix->row_ptr[i - 1];
    }

    aMatrix->row_ptr[0] = 0;

    // free the temporary array
    free(elements);
}

void CreateNonConstantMatrix(const CSRMatrix aMatrix, CSRMatrix *nonConsantMatrix)
{
    nonConsantMatrix->num_rows = aMatrix.num_rows;
    nonConsantMatrix->num_cols = aMatrix.num_cols;
    nonConsantMatrix->num_non_zeros = aMatrix.num_non_zeros;
    // allocate memory for the matrix
    nonConsantMatrix->csr_data = (double *)malloc(nonConsantMatrix->num_non_zeros * sizeof(double));
    nonConsantMatrix->col_ind = (int *)malloc(nonConsantMatrix->num_non_zeros * sizeof(int));
    nonConsantMatrix->row_ptr = (int *)malloc((nonConsantMatrix->num_rows + 1) * sizeof(int));
    
    if (nonConsantMatrix->csr_data == NULL || nonConsantMatrix->col_ind == NULL || nonConsantMatrix->row_ptr == NULL)
    {
        printf("Error: memory allocation failed\n");
        exit(1);
    }

    memcpy(nonConsantMatrix->csr_data, aMatrix.csr_data, nonConsantMatrix->num_non_zeros * sizeof(double));
    memcpy(nonConsantMatrix->col_ind, aMatrix.col_ind, nonConsantMatrix->num_non_zeros * sizeof(int));
    memcpy(nonConsantMatrix->row_ptr, aMatrix.row_ptr, (nonConsantMatrix->num_rows + 1) * sizeof(int));
    // we need to copy as if we don't then the original matrix will be changed

    CSRTranspose(nonConsantMatrix);
    printf("Non constant matrix created\n");
}
