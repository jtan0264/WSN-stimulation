#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#include "structure.h"
#include "calculation.h"

#define PI 3.14159265358979323846

/****************************************************************************************
computes the distance between the two points of given latitudes and longitudes
the function takes in 2 pairs of latitude and longitude values 
this function has been obtained from GeoDataStructure
source: https://www.geodatasource.com/developers/c
*****************************************************************************************/
float distance(float lat1, float lon1, float lat2, float lon2, char unit) {
  float theta, dist;
  if ((lat1 == lat2) && (lon1 == lon2)) {
    return 0;
  }
  else {
    theta = lon1 - lon2;
    dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
    dist = acos(dist);
    dist = rad2deg(dist);
    dist = dist * 60 * 1.1515;
    switch(unit) {
      case 'M':
        break;
      case 'K':
        dist = dist * 1.609344;
        break;
      case 'N':
        dist = dist * 0.8684;
        break;
    }
    return (dist);
  }
}

// this funciton converts decimal degrees to radians
float deg2rad(float deg) {
  return (deg * PI / 180);
}

// this function converts radians to decimal degrees
float rad2deg(float rad) {
  return (rad * 180 / PI);
}

// computes whether the reading of a specific node constitues for an earthquake when compared with its neighbourign nodes
bool compare_readings(reading_struct reading_a, reading_struct reading_b, double dist_diff_thresh, double mag_diff_thresh) {
  float distance_between = distance(reading_a.latitude, reading_a.longitude, reading_b.latitude, reading_b.longitude, 'K');
  float mag_between = fabsf(reading_a.magnitude - reading_b.magnitude);

  if (distance_between < dist_diff_thresh && mag_between < mag_diff_thresh) {
    return true;    // passes as a possible earthquake
  }
  else {
    return false;   // reading doesn't match with its neighbours
  }
}