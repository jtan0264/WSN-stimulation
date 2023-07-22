#include <mpi.h>

void run_seafloor_nodes(MPI_Comm world_comm, MPI_Comm sea_comm, int dims[]);

void* request_neighbour_reading(void *arg);

void* send_reading(void *arg);