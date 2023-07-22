extern double mag_thresh;
extern double dist_diff_thresh;	// let the distance be in the unit kilometers
extern double mag_diff_thresh;
extern double balloon_thresh;

typedef struct {
    int year;
    int month; 
    int day;
    int hour;
    int minute;
    int second; 
    float latitude;
    float longitude;
    float magnitude;
    float depth;
    int rank;
    int x_coor;
    int y_coor;
    int actual_neighbours;
} reading_struct;

typedef struct {
  int sentinel_value;
  double threshold;
} balloon_struct;

reading_struct node_reading(int m, int n, int coord[], int my_rank);

float generate_random_float(char type, int m, int n, int coord[], int my_rank);

reading_struct balloon_reading(double balloon_thresh);

float generate_random_float_balloon(char type, int threshold);