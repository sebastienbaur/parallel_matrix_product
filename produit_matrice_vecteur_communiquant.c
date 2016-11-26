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

void fill(double *vector_to_be_filled, double *filling_vector){
    int size = sizeof(filling_vector);
    int i;
    for (i=0;i<size;i++){
        vector_to_be_filled[i] = filling_vector[i];
    }
}

int main(int argc, char *argv[]) {
    //Start MPI...
    MPI_Init(&argc, &argv);
    int rank; /* rank of the process */
    int size; /* number of processes */
    int source; /* rank of the sender */
    int dest; /* rank of the receiver */
    int tag;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // get data from the command line
    int column_number = atoi(argv[1]);
    int line_number = atoi(argv[2]);

    // the first processor won't be used to compute the product, so that you have size - 1 processors left
    int slice_number = size - 1;
    //Initialize local result and matrices
    double *local_result = NULL;
    double *local_matrix = NULL;
    double *local_vector=NULL;
    int slice_size;
    int max_received_size = (line_number%slice_number) + line_number/slice_number; // process 0 receives subvector of size at most max_received_size


    /* processor 0 generates a matrix and a vector,
    slices each and sends slices to the appropriate processors*/
    if(rank==0){
        slice_size = line_number/slice_number;
        int min_value_in_matrix = atoi(argv[3]);
        int max_value_in_matrix = atoi(argv[4]);
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
            local_vector = vector + (dest-1)*slice_size;
            MPI_Send(local_vector, slice_size, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
            MPI_Send(&slice_size, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
            MPI_Send(local_matrix, slice_size*column_number, MPI_DOUBLE, dest, 3, MPI_COMM_WORLD);
        }
        // the last processor should be treated differently since it doesn't always have the same number of rows
        dest = size-1;
        local_matrix = matrix +  column_number*(dest-1)*slice_size;
        local_vector = vector + (dest-1)*slice_size;
        slice_size = (line_number%slice_number) + line_number/slice_number; /* not very smart, especially if line_number%slice_number = line_number - 1. But it's a beginning :)*/
        MPI_Send(local_vector, slice_size, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
        MPI_Send(&slice_size, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
        MPI_Send(local_matrix, slice_size*column_number, MPI_DOUBLE, dest, 3, MPI_COMM_WORLD);
        free(matrix);
        free(vector);


        // Transmits data from a processor to each of the other processors
        int received_size;
        double *received_vector;
        received_vector = malloc(sizeof(double)*max_received_size);
        for(source=1;source<size;source++){
            MPI_Recv(&received_size, 1, MPI_INT, source,1,MPI_COMM_WORLD,&status);
            MPI_Recv(received_vector, max_received_size, MPI_DOUBLE, source, 2, MPI_COMM_WORLD,&status);
            //printf("Vector received from processor %d \n", source);
            //print_matrix(received_vector, received_size, 1);
            //printf("\n");
            for(dest=1;dest<size;dest++){
                if(dest!=source){
                    MPI_Send(&received_size, 1, MPI_INT, dest, 2*source, MPI_COMM_WORLD);
                    MPI_Send(received_vector, max_received_size, MPI_DOUBLE, dest, 2*source+1, MPI_COMM_WORLD);
                }
            }
        }
        free(received_vector);


        // Get all the subresults and prints the global result
        double *global_result;
        global_result = malloc(sizeof(double)*line_number);
        local_result = malloc(sizeof(double)*max_received_size);
        source = 1;
        int compteur;
        for(source=1;source<size;source++){
            MPI_Recv(&slice_size, 1, MPI_INT, source,1,MPI_COMM_WORLD,&status);  // received before the local result
            MPI_Recv(local_result, slice_size, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &status);
            printf("Received a vector of size %d from processor %d \n",slice_size,source);
            print_matrix(local_result, slice_size, 1);
            for(compteur=0; compteur<slice_size;compteur++){
                global_result[(source-1)*(line_number/slice_number)+compteur] = local_result[compteur];
                printf("\n Global_result at index %d : %f\n", (source-1)*slice_size+compteur,local_result[compteur]);
            }
        }
        printf("\n Result : \n");
        print_matrix(global_result, line_number, 1);
        free(global_result);
        free(local_result);
    }


    else{ // rank != 0
        /* Receives all sliced data from processor 0 :
        local vector and local matrix, as well as slice_size*/
        source=0;
        dest=0;
        MPI_Recv(&slice_size, 1, MPI_INT, source,1,MPI_COMM_WORLD,&status); // receive first the size (i.e the number of rows of the considered slice)
        local_vector = malloc(slice_size*sizeof(double));
        MPI_Recv(local_vector, slice_size, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &status);
        local_matrix = malloc(sizeof(double)*column_number*slice_size);
        MPI_Recv(local_matrix, column_number*slice_size, MPI_DOUBLE, source, 3, MPI_COMM_WORLD, &status);

        // Initialises local_result
        local_result = malloc(sizeof(double)*slice_size);
        int k;
        for(k=0;k<slice_size;k++){
            local_result[k] = 0;
        }


        // Creating data that will be exchanged between the processors
        double *received_vector;
        int received_size = slice_size;
        received_vector = malloc(sizeof(double)*max_received_size);
        fill(received_vector, local_vector);
        // Sends local_vector to processor 0 that will transmit it to the other processors
        MPI_Send(&received_size, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
        MPI_Send(received_vector, max_received_size, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD); // sent before the local result
        // Receives the other local vectors and uses them for the computation
        int j;
        int i;
        for(tag=1;tag<size;tag++){ // receives the local vectors from processor 0 and computes the subproducts
            if(tag!=rank){
                MPI_Recv(&received_size, 1, MPI_INT, source, 2*tag, MPI_COMM_WORLD, &status);
                MPI_Recv(received_vector, max_received_size, MPI_DOUBLE, source, 2*tag+1, MPI_COMM_WORLD, &status);
                for(i=0;i<slice_size;i++){
                    for(j=0;j<received_size;j++){
                        local_result[i] += local_matrix[i*column_number + (tag-1)*slice_size + j]*received_vector[j];
                    }
                }
            }
        }
        for(i=0;i<slice_size;i++){ // Computes the subproduct for its own local vector
            for(j=0;j<slice_size;j++){
                local_result[i] += local_matrix[i*column_number + (rank-1)*slice_size + j]*local_vector[j];
            }
        }


        // Sends the local_result to the processor 0 so that it will merge them
        MPI_Send(&slice_size, 1, MPI_INT, dest, 1, MPI_COMM_WORLD); // sent before the local result
        MPI_Send(local_result, slice_size, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
        free(local_matrix);
        free(local_result);
        free(local_vector);
        free(received_vector);
    }
    MPI_Finalize();
    return 0;
}
