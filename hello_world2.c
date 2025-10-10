#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int message = 0;

    if (world_rank == 0) {
        printf("Hello World from rank %d out of %d message %d \n", world_rank, world_size, message);
        if (world_size > 1) {
            MPI_Send(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(&message, 1, MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        message = world_rank;
        printf("Hello World from rank %d out of %d message %d \n", world_rank, world_size, message);
        if (world_rank < world_size - 1) {
            MPI_Send(&message, 1, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    return 0;
}
