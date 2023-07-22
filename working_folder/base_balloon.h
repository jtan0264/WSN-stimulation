#include <mpi.h>
#include<stdbool.h>

void run_base_station(MPI_Comm, MPI_Comm, double, double, double);

void *balloon_sensor(void *pArg);

void *send_sentinel();

void *receive_reading();

void write_log_conclusive(FILE *logFile, int iteration, char conclusion[], reading_struct node_report, reading_struct balloon_report, reading_struct neighbour_reading[]);

void write_log_inconclusive(FILE *logFile, int iteration, char conclusion[], reading_struct node_report, reading_struct neighbour_reading[]);

void write_summary(FILE *logFile, int events, int conclusive_count, int inconclusive_count);