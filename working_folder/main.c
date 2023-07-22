#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include<stdbool.h>

#include "structure.h"
#include "calculation.h"
#include "seafloor_nodes.h"
#include "base_balloon.h"

// global variables
MPI_Datatype TypeReading;	// allows other files to use this datatype
double mag_thresh;
double dist_diff_thresh;	// let the distance be in the unit kilometers
double mag_diff_thresh;
double balloon_thresh;

int main(int argc, char **argv) {
	// variables for the main MPI communicator MPI_COMM_WORLD
	int my_rank, size, provided;

	// variable for storing the dimension of seafloor nodes cartesian
	int dims[2];
	int m, n;   	// m x n nodes in a cartesian grid
	char *eptr;		// variable for the argument of strtod function

	// communicator variable for base station (master) and seafloor nodes (slaves)
	MPI_Comm new_comm;

    // start up initial MPI environment
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /************************************************
    Initialising the reading_struct in the main file
    ************************************************/
    // variables for the reading struct
    int num_attr = 14;
    MPI_Datatype type[14] ={MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    int blocklen[14] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint disp[14];

    // storing the addresses of the struct variables into MPI_Aint disp[10] array
    disp[0] = offsetof(reading_struct, year);
    disp[1] = offsetof(reading_struct, month);
    disp[2] = offsetof(reading_struct, day);
    disp[3] = offsetof(reading_struct, hour);
    disp[4] = offsetof(reading_struct, minute);
    disp[5] = offsetof(reading_struct, second);
    disp[6] = offsetof(reading_struct, latitude);
    disp[7] = offsetof(reading_struct, longitude);
    disp[8] = offsetof(reading_struct, magnitude);
    disp[9] = offsetof(reading_struct, depth);
    disp[10] = offsetof(reading_struct, rank);
    disp[11] = offsetof(reading_struct, x_coor);
    disp[12] = offsetof(reading_struct, y_coor);
	disp[13] = offsetof(reading_struct, actual_neighbours);

    MPI_Type_create_struct(num_attr, blocklen, disp, type, &TypeReading);
    MPI_Type_commit(&TypeReading);

    /*******************************************************************************************************************
	m and n values are taken as command line arguments
	the threshold values of magnitude, distance difference and magnitude difference are also taken as command line input
	the input is in the order
		m	n	mag_thresh	dist_diff_thresh	mag_diff_thresh		balloon_thresh
	********************************************************************************************************************/
	if (argc == 7) {
		m = atoi(argv[1]);
		n = atoi(argv[2]);
		mag_thresh = strtod(argv[3], &eptr);
		dist_diff_thresh = strtod(argv[4], &eptr);
		mag_diff_thresh = strtod(argv[5], &eptr);
		balloon_thresh = strtod(argv[6], &eptr);

		dims[0] = m;
		dims[1] = n;

		if ((m * n) != (size - 1)) {
			if (my_rank == 0) {
				printf("There are not enough processes to create the number of nodes\n");
			}
			MPI_Finalize();
			return 0;
		}
	}

	else {	
		if (my_rank == 0) {
			printf("User has to input values  1)m\t2)n\t3)mag_thresh\t4)dist_diff_thresh\t5)mag_diff_thresh\t6)balloon_thresh (total of 6 values)\n");
		}
		MPI_Finalize();
		return 0;
	}

	// splitting the base station process and seafloor nodes processes
	MPI_Comm_split(MPI_COMM_WORLD, my_rank == 0, 0, &new_comm); // color will either be 0 or 1 

	/*******************************************************
	this would be the function call to initiate the program
	*******************************************************/
	if (my_rank == 0) {
		// call base_station function (JiaSheng)
		run_base_station(MPI_COMM_WORLD, new_comm, dist_diff_thresh, mag_diff_thresh, balloon_thresh);
		printf("EXIT RUN_BASE_STATION\n");
	}
	else {
		// call seafloor_nodes function (KangZhuang)
		run_seafloor_nodes(MPI_COMM_WORLD, new_comm, dims);
		printf("EXIT RUN_SEAFLOOR_NODES\n");
	}
    
    /**********************************************************
	cleaning up the memory
	**********************************************************/
	// clean up the type 
	MPI_Type_free(&TypeReading);

	MPI_Finalize();
	return 0;
}