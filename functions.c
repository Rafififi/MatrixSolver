#include <stdio.h>
#include <stdlib.h>
#include "functions.h"


void ReadMMtoCSR(const char *filename, CSRMatrix *aMatrix)
{
    FILE *fileName = fopen(filename, "r");
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
    
    int *rows;
    int *cols;
    double *data;

    rows = (int *)calloc(aMatrix->num_non_zeros, sizeof(int));
    cols = (int *)calloc(aMatrix->num_non_zeros,  sizeof(int));
    data = (double *)calloc(aMatrix->num_non_zeros,  sizeof(double));
    for (int i = 0; i<aMatrix->num_non_zeros; i++)
    {
        fscanf(fileName, "%d %d %lf", &rows[i], &cols[i], &data[i]);
        c = getc(fileName);
    }
    printf("\n");


    aMatrix->row_ptr = (int *)calloc(aMatrix->num_rows, sizeof(int));
    aMatrix->csr_data = (double *)calloc(aMatrix->num_non_zeros, sizeof(double));
    aMatrix->col_ind = (int *)calloc(aMatrix->num_cols, sizeof(int));

    //sort the rows and apply the same permutation to cols and data then sort cols making sure that the data and rows are sorted with it and then cols are sorted with the rows and data

    //sort the rows
    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        for (int j = i+1; j < aMatrix->num_non_zeros; j++)
        {
            if (rows[i] > rows[j])
            {
                int temp = rows[i];
                rows[i] = rows[j];
                rows[j] = temp;

                int temp2 = cols[i];
                cols[i] = cols[j];
                cols[j] = temp2;

                double temp3 = data[i];
                data[i] = data[j];
                data[j] = temp3;
            }

            if (rows[i] == rows[j])
            {
                if (cols[i] > cols[j])
                {
                    int temp = cols[i];
                    cols[i] = cols[j];
                    cols[j] = temp;

                    double temp2 = data[i];
                    data[i] = data[j];
                    data[j] = temp2;
                }
            }
        }
    }    

    //count how many values there are in each row and store it in row_ptr
    // all the values in row_ptr need to be equal to that number and the numbers before it added together
    int count = 0;
    int row = 0;
    for (int i = 0; i <= aMatrix->num_non_zeros; i++)
    {
        if (rows[i] == row)
        {
            count++;

        }
        else
        {
            aMatrix->row_ptr[row] = count;
            count++;
            row++;
        }
        double temp1 = data[i];
        aMatrix->csr_data[i] = temp1;

        int temp2 = cols[i];
        aMatrix->col_ind[i] = temp2 - 1;
    }    

    printf("Number of non-zero elements: %d\n", aMatrix->num_non_zeros);
    printf("Row Pointer:");
    for (int i = 0; i <= aMatrix->num_rows; i++)
    {
        printf(" %d", aMatrix->row_ptr[i]);
    }

    printf("\n");
    printf("Coloumn Index:");
    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        printf(" %d", aMatrix->col_ind[i]);
    }
    printf("\n");
    printf("Values:");
    for (int i = 0; i < aMatrix->num_non_zeros; i++)
    {
        printf(" %lf", aMatrix->csr_data[i]);
    }
    printf("\n");
    
    free(rows);
    free(cols);
    free(data);
    fclose(fileName);

}

void spmv_csr(const CSRMatrix *AMatrix, const double *x, double *y)
{
    double *product = (double *)calloc(AMatrix->num_non_zeros, sizeof(double));
    for (int row = 0; row < AMatrix->num_rows; row++)
    {

        for (int j = AMatrix->row_ptr[row]; j < AMatrix->row_ptr[row + 1]; j++)
        {
            product[row] += AMatrix->csr_data[j] * x[AMatrix->col_ind[j]];
        }
        y[row] = product[row];
    }
    printf("Product:");
    for (int i = 0; i < AMatrix->num_rows; i++)
    {
        printf(" %lf", y[i]);
    }
}



