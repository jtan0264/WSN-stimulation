#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include<stdbool.h>

#include "structure.h"
#include "base_balloon.h"
#include "calculation.h"
#include "seafloor_nodes.h"

// symbolic constants
#define SENTINEL 1
#define CYCLE 2
#define BALLOON_ARRAY_ROW 10

// symbolic constants of MPI tag values 
#define REQUEST_TAG 1
#define RETURN_TAG 2
#define SENTINEL_TAG 3
#define REPORT_TAG 4
#define NEIGHBOUR_TAG 5

// global variables
MPI_Datatype TypeReading;
pthread_mutex_t g_Mutex = PTHREAD_MUTEX_INITIALIZER;
reading_struct *report_global_array = NULL;     // global array of readings for the balloon

int iteration = 0;
int conclusive_count = 0;
int inconclusive_count = 0;

// create a log file
FILE *logFile;

// mutex 
pthread_mutex_t global_array_mutex = PTHREAD_MUTEX_INITIALIZER;

void run_base_station(MPI_Comm world_comm, MPI_Comm comm, double dist_diff_thresh, double mag_diff_thresh, double balloon_thresh) {
    // initializing the global array for reading
    report_global_array = (reading_struct *)malloc(BALLOON_ARRAY_ROW * sizeof(reading_struct));

    // // initializing all values of the global array 
    reading_struct reading;
    reading.year = 0;
    reading.month = 0;
    reading.day = 0;
    reading.hour = 0;
    reading.minute = 0;
    reading.second = 0;
    reading.latitude = -999999;
    reading.longitude =-999999;
    reading.magnitude = -999999;
    reading.depth = -999999;
    for(int i = 0; i< BALLOON_ARRAY_ROW; i++){
        report_global_array[i] = reading;
    }

    logFile = fopen("log.txt", "w+"); 
    // test wheter log file exists
    if (logFile == NULL) {
        printf("Error! Could not open the sentinel text file\n");
        exit(-1); // must include stdlib.h
    }

    // values that will be passed as thread function argument
    // store the value dist_diff_thresh and mag_diff_thresh into differences array as we would be passing this array as a thread function argument
    double *differences = (double*)malloc((2)*sizeof(double));
    differences[0] = dist_diff_thresh;
    differences[1] = mag_diff_thresh;

    balloon_struct b_arg_package;
    b_arg_package.sentinel_value = 0;
    b_arg_package.threshold = balloon_thresh;

    // thread id
    pthread_t tid0, tid1, tid2;

    // intiialising the mutex
    pthread_mutex_init(&global_array_mutex, NULL);

    // create the thread to run the thread functions
    pthread_create(&tid0, 0, balloon_sensor, &b_arg_package); 
    pthread_create(&tid2, NULL, receive_reading, (void*)differences); 

    int sentinel_value = 0;

    do {
        FILE *in_file = fopen("termination.txt", "r"); // read only

        // test for files not existing.
        if (in_file == NULL) {
            printf("Error! Could not open the sentinel text file\n");
            exit(-1); // must include stdlib.h
        }

        // read from file/keyboard. remember the ampersands!
        fscanf(in_file, "%d", &sentinel_value);

        if (sentinel_value == SENTINEL) {
            pthread_cancel(tid2);
            b_arg_package.sentinel_value = 1;

            // pthread_join(tid0, NULL);
            pthread_cancel(tid0);

            pthread_mutex_destroy(&global_array_mutex);
            fclose(in_file);
            pthread_create(&tid1, NULL, send_sentinel, NULL); // thread function to send the termination sentinel value to the seafloor nodes
        }

        pthread_join(tid1, NULL);
        sleep(CYCLE);
    } while (sentinel_value != SENTINEL);

    // write the summary
    write_summary(logFile, iteration, conclusive_count, inconclusive_count);

    // close the log file
    fclose(logFile);
    
}

/******************************************************
the threaed function that act as the balloon sensor
******************************************************/
void *balloon_sensor(void *pArg) {
    reading_struct reading;
    int pointer = 0;
    balloon_struct *paramPtr = (balloon_struct *)pArg;
    
    // unpacking the argument pointer
    int sentinel = paramPtr->sentinel_value;
    int threshold = paramPtr->threshold;

    while (sentinel != SENTINEL) {
        sentinel = paramPtr->sentinel_value;

        // generate balloon reading
        reading = balloon_reading(threshold);

        // placing the reading into the queue
        pthread_mutex_lock(&global_array_mutex);
        report_global_array[pointer] = reading;
        pthread_mutex_unlock(&global_array_mutex);

        // this allows for FIFO
        pointer++;
        pointer = pointer % BALLOON_ARRAY_ROW;
        
        sleep(CYCLE);
    }
    printf("Balloon Sensor Thread terminates\n");
    fflush(stdout);

    return 0;
}

/******************************************************************************
thread function that sends out the termination message to the seafloor nodes
******************************************************************************/
void *send_sentinel() {
    int sentinel_value = 1;
    int world_size;     // world_size are all the MPI processes in this program (base station + seafloor nodes)
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Request request;

    // the seafloor nodes are from Rank 1 to Rank world_size -1
    for (int i = 1; i < world_size; i++) {
        MPI_Isend(&sentinel_value, 1, MPI_INT, i, SENTINEL_TAG, MPI_COMM_WORLD, &request);
    }
    pthread_exit(NULL);
}

/*****************************************************************************
 thread function that handles the receiving of report from seafloor nodes
 this would free up the base station to carry out other tasks
 *****************************************************************************/
void *receive_reading() {
    reading_struct node_report;
    reading_struct balloon_report;
    int running = 1;
    bool match = false;
    bool event_conclusive = false;
    // int iteration = 0;

    MPI_Request request;
    MPI_Status status;

    while (running) { 
        event_conclusive = false;

        // receiving the reading from the seafloor nodes
        MPI_Recv(&node_report, 1, TypeReading, MPI_ANY_SOURCE, REPORT_TAG, MPI_COMM_WORLD, &status);
        printf("(BASE STATION)\t REPORT SENT FROM SEANODE %d\t Magnitude: %lf\n", status.MPI_SOURCE -1, node_report.magnitude);
        fflush(stdout);

        // receiving the neighbour reading of the seafloor_nodes
        reading_struct neighbour_reading[node_report.actual_neighbours];
        for (int i = 0; i < node_report.actual_neighbours; i++) {
            MPI_Recv(&neighbour_reading[i], 1, TypeReading, status.MPI_SOURCE, NEIGHBOUR_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // NEIGHBOUR


        // compares the received reading with the global array of readings
        for (int i = 0; i < BALLOON_ARRAY_ROW; i++) {
            match = compare_readings(node_report, report_global_array[i], dist_diff_thresh, mag_diff_thresh);
            
            if (match) {   // found that the 2 readings match, hence we would log this event as conclusive
                event_conclusive = true; 
                balloon_report = report_global_array[i];
                break;
            }
        }

        // log iteration
        iteration++;

        // a conclusive event
        if (event_conclusive) {
            char conclusive[15] = "CONCLUSIVE";
            conclusive_count++;
            // log to file
            write_log_conclusive(logFile, iteration, conclusive, node_report, balloon_report, neighbour_reading);
        }
        // an inconclusive event 
        else { 
            char inconclusive[15] = "INCONCLUSIVE";
            inconclusive_count++;
            // log to file
            write_log_inconclusive(logFile, iteration, inconclusive, node_report, neighbour_reading);
        }
    }
}

/****************************************************
function that writes conclusive event to log file
*****************************************************/
void write_log_conclusive(FILE *logFile, int iteration, char conclusion[], reading_struct node_report, reading_struct balloon_report, reading_struct neighbour_reading[]) {
    // get current time for logged time
    struct tm *sTm;
    time_t now = time(0);
    sTm = localtime(&now);
    
    fprintf(logFile, "--------------------------------------------------------------------------------------------\n");
    fprintf(logFile, "%-10s %d\n", "Iteration:", iteration);
    fprintf(logFile, "%-32s %d/%02d/%02d\t %02d:%02d:%02d\n", "Logged Time:", sTm -> tm_year + 1900, sTm -> tm_mon, sTm -> tm_mday, sTm-> tm_hour, sTm-> tm_min, sTm-> tm_sec);
    fprintf(logFile, "%-32s %d/%02d/%02d\t %02d:%02d:%02d\n", "Alert Reported Time:", node_report.year, node_report.month, node_report.day, node_report.hour, node_report.minute, node_report.second);
    fprintf(logFile, "%-10s %s\n", "Alert type:", conclusion);

    fprintf(logFile, "\n%s %35s %27s\n", "Reporting Node", "Seismic Coordinate", "Magnitude");
    fprintf(logFile, "%-4d%2s(%d, %d)%23s(%.2f, %.2f) %20s %.2f\n", node_report.rank, "", node_report.x_coor, node_report.y_coor, "", node_report.latitude, node_report.longitude, "", node_report.magnitude);

    for (int i = 0; i < node_report.actual_neighbours; i++) {
        fprintf(logFile, "\nAdjacent Node\n");
        fprintf(logFile, "%d\t (%d, %d)\n", neighbour_reading[i].rank, neighbour_reading[i].x_coor,neighbour_reading[i].y_coor);
        fprintf(logFile, "%50s  %s  %s   %s\n","Seismic Coordinate", "Coord Diff(km)", "Magnitude", "Magnitude Diff");
        fprintf(logFile, "%33s%.2f, %.2f%s %14.2f %14.2f %14.2f\n", "(", neighbour_reading[i].latitude, neighbour_reading[i].longitude, ")", distance(node_report.latitude, node_report.longitude, neighbour_reading[i].latitude, neighbour_reading[i].longitude, 'K'), neighbour_reading[i].magnitude, fabs(node_report.magnitude - neighbour_reading[i].magnitude));
    }
    
    fprintf(logFile, "\n%-32s %d/%02d/%02d\t %02d:%02d:%02d\n", "Balloon Seismic Reported Time:", balloon_report.year, balloon_report.month, balloon_report.day, balloon_report.hour, balloon_report.minute, balloon_report.second);
    fprintf(logFile, "%s (%.2f, %.2f)\n", "Balloon Seismic Reported Coord:", balloon_report.latitude, balloon_report.longitude);
    fprintf(logFile, "%s %.2f\n", "Coord Diff between Balloon and Reporting Node:", distance(node_report.latitude, node_report.longitude, balloon_report.latitude, balloon_report.longitude, 'K'));
    fprintf(logFile, "%s %.2f\n", "Magnitude Diff between Balloon and Reporting Node:", fabsf(node_report.magnitude - balloon_report.magnitude));
    fprintf(logFile, "--------------------------------------------------------------------------------------------\n\n");
    
}

/***************************************************
function that writes inconclusive event to log file
***************************************************/
void write_log_inconclusive(FILE *logFile, int iteration, char conclusion[], reading_struct node_report, reading_struct neighbour_reading[]) {
    // get current time for logged time
    struct tm *sTm;
    time_t now = time(0);
    sTm = localtime(&now);

    fprintf(logFile, "--------------------------------------------------------------------------------------------\n");
    fprintf(logFile, "%-10s %d\n", "Iteration:", iteration);
    fprintf(logFile, "%-32s %d/%02d/%02d\t %02d:%02d:%02d\n", "Logged Time:", sTm -> tm_year + 1900, sTm -> tm_mon, sTm -> tm_mday, sTm-> tm_hour, sTm-> tm_min, sTm-> tm_sec);
    fprintf(logFile, "%-32s %d/%02d/%02d\t %02d:%02d:%02d\n", "Alert Reported Time:", node_report.year, node_report.month, node_report.day, node_report.hour, node_report.minute, node_report.second);
    fprintf(logFile, "%-10s %s\n", "Alert type:", conclusion);

    fprintf(logFile, "\n%s %35s %27s\n", "Reporting Node", "Seismic Coordinate", "Magnitude");
    fprintf(logFile, "%-4d%2s(%d, %d)%23s(%.2f, %.2f) %20s %.2f\n", node_report.rank, "", node_report.x_coor, node_report.y_coor, "", node_report.latitude, node_report.longitude, "", node_report.magnitude);
    
    for (int i = 0; i < node_report.actual_neighbours; i++) {
        fprintf(logFile, "\nAdjacent Node\n");
        fprintf(logFile, "%d\t (%d, %d)\n", neighbour_reading[i].rank, neighbour_reading[i].x_coor,neighbour_reading[i].y_coor);
        fprintf(logFile, "%50s  %s  %s   %s\n","Seismic Coordinate", "Coord Diff(km)", "Magnitude", "Magnitude Diff");
        fprintf(logFile, "%33s%.2f, %.2f%s %14.2f %14.2f %14.2f\n", "(", neighbour_reading[i].latitude, neighbour_reading[i].longitude, ")", distance(node_report.latitude, node_report.longitude, neighbour_reading[i].latitude, neighbour_reading[i].longitude, 'K'), neighbour_reading[i].magnitude, fabs(node_report.magnitude - neighbour_reading[i].magnitude));
    }

    fprintf(logFile, "\nNo Balloon Seismic Reading as the earthquake event is inconclusive\n");
    fprintf(logFile, "--------------------------------------------------------------------------------------------\n\n");
}

void write_summary(FILE *logFile, int events, int conclusive_count, int inconclusive_count) {
    fprintf(logFile,"Summary of earthquakes\n");
    fprintf(logFile,"Number of Events: %d\n", iteration);
    fprintf(logFile,"Number of Conclusive Events: %d\n", conclusive_count);
    fprintf(logFile,"Number of Inconclusive Events: %d\n", inconclusive_count);
}


