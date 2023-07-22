#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>
#include <stdio.h>
#include <pthread.h>
#include <stdbool.h>

#include "seafloor_nodes.h"
#include "structure.h" 
#include "calculation.h"

// symbolic constants
#define SENTINEL 1
#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1
#define NUM_NEIGHBOURS 4
#define CYCLE 2

// symbolic constants of MPI tag values 
#define REQUEST_TAG 1
#define RETURN_TAG 2
#define SENTINEL_TAG 3
#define REPORT_TAG 4
#define NEIGHBOUR_TAG 5

// global variable
MPI_Comm comm2D;	// the communicator that would hold the new cartersian grid
pthread_mutex_t request_mutex = PTHREAD_MUTEX_INITIALIZER;	// mutex to give threads the ability to lock 
MPI_Datatype TypeReading;
reading_struct reading;
reading_struct recv_buf[NUM_NEIGHBOURS];
int actual_neighbours = 0;
int cycle_count = 0;

/*
the function that is called by main.c, this function would simulate the behaviour of seanodes
arguments
    world_comm - communicator including the base station
    sea_comm - communicator that involves the seafloor nodes only
*/
void run_seafloor_nodes(MPI_Comm world_comm, MPI_Comm sea_comm, int dims[]) {
    // cartesian topology variables
	int ndims=2, size, my_rank, global_rank, reorder, my_cart_rank, ierr, sentinel_value;
	int top_nbr, bottom_nbr;
	int left_nbr, right_nbr;
	int *neighbours = NULL;	
	int coord[ndims];
	int wrap_around[ndims];
	char *eptr;		// variable for the argument of strtod function

	// threads id
	pthread_t request_thread;
	pthread_t send_reading_thread;
	
    // getting the attributes of the seafloor nodes communicator
    MPI_Comm_size(sea_comm, &size);
    MPI_Comm_rank(sea_comm, &my_rank); 
	MPI_Comm_rank(world_comm, &global_rank); 

    /*********************************************
	create cartesian topology for processes 
    *********************************************/
	MPI_Dims_create(size, ndims, dims);

	/*************************************************
    create cartesian mapping 
    *************************************************/
	wrap_around[0] = 0;
	wrap_around[1] = 0;     // periodic shift is false in both dimensions
	reorder = 1;
	ierr =0;
	ierr = MPI_Cart_create(sea_comm, ndims, dims, wrap_around, reorder, &comm2D);
	if(ierr != 0) {
        printf("ERROR[%d] creating CART\n",ierr);
    }
	
	// find my coordinates in the cartesian communicator group 
	MPI_Cart_coords(comm2D, my_rank, ndims, coord); // coordinated is returned into the coord array

	// use my cartesian coordinates to find my rank in cartesian group
	MPI_Cart_rank(comm2D, coord, &my_cart_rank);    // the opposite function of MPI_Cart_coords()

	// creates the communicator that allows nodes to only communicate with its adjacent nodes
	MPI_Cart_shift( comm2D, SHIFT_ROW, DISP, &top_nbr, &bottom_nbr );
	MPI_Cart_shift( comm2D, SHIFT_COL, DISP, &left_nbr, &right_nbr );

	// put the neighbouring ranks into an array
	int num_neighbours = 4;
	// int neighbours[num_neighbours];
	neighbours = (int*) malloc (num_neighbours * sizeof(int));	// array to store the neighbourign process rank 

	// non existing neighbours would be set as value -2
	neighbours[0] = top_nbr;
	neighbours[1] = bottom_nbr;
	neighbours[2] = left_nbr;
	neighbours[3] = right_nbr;

	// gettting the actual number of neighbours as not all nodes have 4 neighbours
	for (int i = 0; i < NUM_NEIGHBOURS; i++) {
		if (neighbours[i] >= 0) {
			neighbours[actual_neighbours] = neighbours[i];
			actual_neighbours ++;
		}
	}	

	/****************************************************************************
	creating thread to handle sending of reading to requesting neighbours
	using threads would free up the process to go on with its other computation
	****************************************************************************/
	// creating mutex 
	pthread_mutex_init(&request_mutex, NULL);
	// thread that handles receiving request from neighbour and sending reading to the requesting neighbour
	pthread_create(&send_reading_thread, 0, send_reading, &comm2D);

    
	// seafloor nodes would continue obtaining reading until the receive the termination message from the base station
	do
	{
		MPI_Request request;
		// check whether the nodes have received the message to terminate		
		MPI_Irecv(&sentinel_value, 1, MPI_INT, 0, SENTINEL_TAG, world_comm, &request);
		int ready;
		MPI_Test(&request, &ready, MPI_STATUS_IGNORE);

		// cycle count variable
		cycle_count ++;
		
    	// reading obtained by nodes is put into the reading_struct
		reading = node_reading(dims[0], dims[1], coord, my_rank);
		reading.actual_neighbours = actual_neighbours;	// set the value of actual_neighbours
		MPI_Barrier(comm2D);	// ensures that all nodes have generated a reading at a particular cycle

		if (reading.magnitude > mag_thresh) {
			printf("(SEAFLOOR NODE) Cycle %d\t Node %d exceeds threshold with magnitude reading of %lf\n", cycle_count, my_rank, reading.magnitude);
			fflush(stdout);
			// array of MPI_Request and MPI_Status 
			MPI_Request send_req[actual_neighbours];
			MPI_Status send_status[actual_neighbours];
			MPI_Request recv_req[actual_neighbours];
			MPI_Status recv_status[actual_neighbours];

			// // array of reading_struct
			// reading_struct recv_buf[actual_neighbours];

			// request_signal sent to neighbours 
			int request_signal = 1;

			for (int i = 0; i < actual_neighbours; i++) {
				MPI_Isend(&request_signal, 1, MPI_INT, neighbours[i], REQUEST_TAG, comm2D, &send_req[i]);
				MPI_Irecv(&recv_buf[i], 1, TypeReading, neighbours[i], RETURN_TAG, comm2D, &recv_req[i]);
			}
			MPI_Waitall(actual_neighbours, send_req, send_status);
			MPI_Waitall(actual_neighbours, recv_req, recv_status);

			// calculate the number of reading matches with neighbours
			int match = 0;
			for (int i = 0; i < actual_neighbours; i++) {
				bool valid = compare_readings(reading, recv_buf[i], dist_diff_thresh, mag_diff_thresh);

				if (valid) {
					match ++;
				}
			}
			
			// a report is sent to the base station 
			if (match >= 2) {
				printf("(SEAFLOOR NODE) Rank %d\t MAGNITUDE: %lf\t REPORT TO BE SENT TO BASE STATION\n", my_rank, reading.magnitude);
				fflush(stdout);
				// blocking send is used here as there is a thread in base_balloon.c that calls the corresponding MPI_Recv
				MPI_Send(&reading, 1, TypeReading, 0, REPORT_TAG, world_comm);
				// another blocking send is used to send all the neighbour reading to the base station
				for (int i = 0; i < actual_neighbours; i++) {
					MPI_Send(&recv_buf[i], 1, TypeReading, 0, NEIGHBOUR_TAG, world_comm);
				}
			}
		}
		sleep(CYCLE);

	} while (sentinel_value != SENTINEL);

	// cleaning up the thread and mutex
	pthread_cancel(send_reading_thread);
	pthread_mutex_destroy(&request_mutex);

	printf("Seafloor Node %d terminates\n", global_rank);
	fflush(stdout);
}

// /*******************************************************************************************************
// this function will be called by a created thread
// where each node would create a thread to handle the receiving of request and sending of its own reading
// using the thread would free up the MPI process to continue with its computation
// argument arg is the pointer to the cartesian communicator
// *******************************************************************************************************/
void* send_reading(void *arg) {
	MPI_Comm comm2D = * (MPI_Comm *)arg;	// obtaining the function argument 

	// local variables
	int probe_flag;
	int request_signal;
	MPI_Request request;
	MPI_Status status;

	while (1) {
		// locking mutex
		pthread_mutex_lock(&request_mutex);

		MPI_Recv(&request_signal, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, comm2D, &status);
		MPI_Send(&reading, 1, TypeReading, status.MPI_SOURCE, RETURN_TAG, comm2D);

		// unlocking the mutex
		pthread_mutex_unlock(&request_mutex);
	}
}