#include <stdbool.h>

float distance(float lat1, float lon1, float lat2, float lon2, char unit);

float deg2rad(float deg);

float rad2deg(float rad);

bool compare_readings(reading_struct reading_a, reading_struct reading_b, double dist_diff_thresh, double mag_diff_thresh);









