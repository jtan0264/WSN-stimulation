#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "structure.h"

/*******************************************************************************************************
simulates the reading of the seafloor nodes by generating random values for the values in the attributes
*******************************************************************************************************/
reading_struct node_reading(int m, int n, int coord[], int my_rank) {
    struct tm *sTm;
    time_t now = time(0);
    sTm = localtime(&now);

    // creating the reading struct and placing the values into the struct
    reading_struct reading;
    reading.year =sTm-> tm_year + 1900;
    reading.month =sTm-> tm_mon;
    reading.day =sTm-> tm_mday;
    reading.hour =sTm-> tm_hour;
    reading.minute =sTm-> tm_min;
    reading.second =sTm-> tm_sec;
    reading.latitude = generate_random_float('L', m, n, coord, my_rank);
    reading.longitude = generate_random_float('l', m, n, coord, my_rank);
    reading.magnitude = generate_random_float('m', m, n, coord, my_rank);
    reading.depth = generate_random_float('d', m, n, coord, my_rank);
	reading.rank = my_rank;
	reading.x_coor = coord[0];
	reading.y_coor = coord[1];

    return reading;
}

/*******************************************************************************************
generates a random number within a certain range given as an input
the function takes in a character value to decide on the range of the randomly generated value
*******************************************************************************************/
float generate_random_float(char type, int m, int n, int coord[], int my_rank) {
    float rand_value;
	float scale; 	// used for generated lat and longi for the cartesian grid
	unsigned int seed = time(NULL) + my_rank;
	int node_range;
	int lower;
	int higher;
	int lat_range = 10;		// these values are set extermely small and not the full range
	int longi_range = 10;

	switch(type) {
		case 'L':	// represents latitude
			node_range = lat_range / m;
			lower = (coord[0] * node_range) - (lat_range/2);
			higher = ((coord[0] + 1) * node_range) - (lat_range/2);
			scale = rand_r(&seed) / (float) RAND_MAX;
			rand_value = lower + scale * (higher - lower);
			break;

		case 'l':	// represents longitude
			node_range = longi_range / n;
			lower = (coord[1] * node_range) - (longi_range/2);
			higher = ((coord[1] + 1) * node_range) - (longi_range/2);
			scale = rand_r(&seed) / (float) RAND_MAX;
			rand_value = lower + scale * (higher - lower);
			break;

		case 'm':	// represents magnitude (greatest ever recorded earthquake: 9.5)
			// rand_value = (float)rand_r(&seed)/(float)(RAND_MAX/10);
			// break;
			higher = 7.0;
			lower = 4.5;
			scale = rand_r(&seed) / (float)RAND_MAX;
			rand_value = lower + scale * (higher - lower);
			break;

		case 'd':	// represents the depth of the earthquake (range 0 to 700km)
			rand_value = (float)rand_r(&seed)/(float)(RAND_MAX/700);
			break;
	}		

	return (rand_value);
}

/*************************************************************************************************************
simulates the reading of the balloon sensor node by generating random values for the values in the attributes
*************************************************************************************************************/
reading_struct balloon_reading(double balloon_thresh) {
	struct tm *sTm;
    time_t now = time(0);
    sTm = localtime(&now);

    // creating the reading struct and placing the values into the struct
    reading_struct reading;
	reading.year =sTm-> tm_year + 1900;
    reading.month =sTm-> tm_mon;
    reading.day =sTm-> tm_mday;
    reading.hour =sTm-> tm_hour;
    reading.minute =sTm-> tm_min;
    reading.second =sTm-> tm_sec;
    reading.latitude = generate_random_float_balloon('L', balloon_thresh);
    reading.longitude = generate_random_float_balloon('l', balloon_thresh);
    reading.magnitude = generate_random_float_balloon('m', balloon_thresh);
    reading.depth = generate_random_float_balloon('d', balloon_thresh);

    return reading;
}

/*******************************************************************************************
generates a random number for the balloon as the balloon has to generate a specific threshold value
the function takes in a character value to decide on the range of the randomly generated value
*******************************************************************************************/
float generate_random_float_balloon(char type, int threshold) {
	float rand_value;
	float scale; // used for generated lat and longi for the cartesian grid
	unsigned int seed = time(NULL);

	int node_range;
	int lower;
	int higher;

	switch (type)
	{
	case 'L': // represents latitude	
		higher = 15;
		rand_value = (float)rand_r(&seed) / (float)(RAND_MAX / higher);
		break;

	case 'l': // represents longitude
		higher = 15;
		rand_value = (float)rand_r(&seed) / (float)(RAND_MAX / higher);
		break;

	case 'm': // represents magnitude (greatest ever recorded earthquake: 9.5)
		higher = 7;
		lower = threshold;
		scale = rand_r(&seed) / (float)RAND_MAX;
		rand_value = lower + scale * (higher - lower);
		break;

	case 'd': // represents the depth of the earthquake (range 0 to 700km)
		higher = 700;
		rand_value = (float)rand_r(&seed) / (float)(RAND_MAX / 700);
		break;
	}

	return (rand_value);
}