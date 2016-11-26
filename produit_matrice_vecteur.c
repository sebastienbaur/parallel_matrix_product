#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>


double random(double min, double max){
    double r = (double)rand()/RAND_MAX;
    return min + r*(max-min);
}

void fill_randomly(double *matrix, int line_number, int column_number, int min_value, int max_value){
    int i = 0;
    int j = 0;
    for(i=0;i<line_number;i++){
        for(j=0;j<column_number;j++){
            //printf("(i,j) = (%d,%d)\n", i, j); /* for debugging purpose */
            //printf("i*column_number + j = %d \n", i*column_number+j); /* for debugging purpose */
            matrix[i*column_number+j] = random(min_value,max_value);
        }
    }
}

void print_matrix(double *matrix, int number_of_rows, int number_of_columns){
    int i;
    int j;
    for(i=0;i<number_of_rows;i++){
        for(j=0;j<number_of_columns;j++){
            printf("%f ", matrix[i*number_of_columns + j]);
        }
        printf("\n");
    }
}

void compute_product(double *matrix, int row_number, int column_number, double *vector, double *result){ /* could be extended to compute matrix matrix product
                                                                                    Doesn't try to manage what happens where the sizes are invalid (too big or too small)...*/
    int i = 0;
    int j = 0;
    for(i=0;i<row_number;i++){
        result[i] = 0;
        for(j=0;j<column_number;j++){
            result[i] += matrix[i*column_number + j]*vector[j];
        }
    }
    printf("\n\n");
}


int main(int argc, char *argv[]) {
    //Start MPI...
    MPI_Init(&argc, &argv);
    int rank; /* rank of the process */
    int size; /* number of processes */
    int source; /* rank of the sender */
    int dest; /* rank of the receiver */
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int column_number = atoi(argv[1]);
    int line_number = atoi(argv[2]);

    // the first processor won't be used to compute the product, so that you have size - 1 processors left
    int slice_number = size - 1;
    //Initialize local result and matrices
    double *local_result = NULL;
    double *local_matrix = NULL;
    int slice_size;


    if(rank==0){
        slice_size = line_number/slice_number;
        int min_value_in_matrix = -12;
        int max_value_in_matrix = 16;
        // generates a matrix
        double *matrix = NULL;
        matrix = malloc(line_number*column_number*sizeof(double));
        fill_randomly(matrix, line_number, column_number, min_value_in_matrix, max_value_in_matrix);
        printf("\n Matrix : \n");
        print_matrix(matrix, line_number, column_number);
        // generates a vector
        double *vector = NULL;
        vector = malloc(line_number*sizeof(double));
        fill_randomly(vector, line_number, 1, min_value_in_matrix, max_value_in_matrix);
        printf("\n Vector : \n");
        print_matrix(vector, line_number, 1);
        for(dest=1;dest<size-1;dest++){ /* the last one will be treated differently for it may have less rows */
            // initializes a local matrix pointer
            local_matrix = matrix +  column_number*(dest-1)*slice_size;
            MPI_Send(vector, line_number, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
            MPI_Send(&slice_size, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
            MPI_Send(local_matrix, slice_size*column_number, MPI_DOUBLE, dest, 3, MPI_COMM_WORLD);
        }
        // the last processor should be treated differently since it doesn't always have the same number of rows
        dest = size-1;
        local_matrix = matrix +  column_number*(dest-1)*slice_size;
        slice_size = (line_number%slice_number) + line_number/slice_number; /* not very smart, especially if line_number%slice_number = line_number - 1. But it's a beginning :)*/
        MPI_Send(vector, line_number, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
        MPI_Send(&slice_size, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
        MPI_Send(local_matrix, slice_size*column_number, MPI_DOUBLE, dest, 3, MPI_COMM_WORLD);
        free(matrix);
        free(vector);


        // Get all the subresults and prints the global result
        double *global_result;
        global_result = malloc(sizeof(double)*line_number);
        int max_received_size = (line_number%slice_number) + line_number/slice_number; // process 0 receives subvector of size at most max_received_size
        local_result = malloc(sizeof(double)*max_received_size);
        source = 1;
        int compteur;
        for(source=1;source<size;source++){
            MPI_Recv(&slice_size, 1, MPI_INT, source,1,MPI_COMM_WORLD,&status);  // received before the local result
            MPI_Recv(local_result, slice_size, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &status);
            for(compteur=0; compteur<slice_size;compteur++){
                global_result[(source-1)*slice_size+compteur] = local_result[compteur];
            }
        }
        printf("\n Result : \n");
        print_matrix(global_result, line_number, 1);
        free(global_result);
        free(local_result);
    }



    else{ // rank != 0
        source=0;
        dest=0;
        double *vector;
        vector = malloc(line_number*sizeof(double));
        MPI_Recv(&slice_size, 1, MPI_INT, source,1,MPI_COMM_WORLD,&status); // receive first the size (i.e the number of rows of the considered slice)
        MPI_Recv(vector, line_number, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &status);
        local_matrix = malloc(sizeof(double)*column_number*slice_size);
        MPI_Recv(local_matrix, column_number*slice_size, MPI_DOUBLE, source, 3, MPI_COMM_WORLD, &status);
        local_result = malloc(sizeof(double)*slice_size);
        compute_product(local_matrix, slice_size, column_number, vector, local_result);
        MPI_Send(&slice_size, 1, MPI_INT, dest, 1, MPI_COMM_WORLD); // sent before the local result
        MPI_Send(local_result, slice_size, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
        free(local_matrix);
        free(local_result);
    }
    MPI_Finalize();
    return 0;
}
