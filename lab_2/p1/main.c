#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "./engine/proj_engine.h"
#include "immintrin.h"


#define TIME(a, b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))
#define TESTN 0
typedef long long ll;

extern int Init(double *data, long long L);
extern int Check(double *data, long long L);

double *cube_blockA;//store data A
double *cube_blockB;//store data B
double *receive_buffer;//store the data received, 6 segments, not all are used
int node_dsize[6]; // this node communicates with 6 nodes at most, this array save the date size to transfer
int node_rank[6]; // save the node rank to communicate
MPI_Request recv_request[6]; // receive request to check needed data has been received. not all are used
MPI_Request send_request[6]; // send requets
int buffer_offset[6]; // the start index of par i to send
int total_nodes;


// record the cube range
typedef struct {
  int offset_x;
  int end_x;
  int offset_y;
  int end_y;
  int offset_z;
  int end_z;
  int nx;
  int ny;
  int nz;
  int x_par;
  int y_par;
  int z_par;
} Info;

int min(int a, int b){
    if (a > b){
        return b;
    } else {
        return a;
    }
}


//Must be modified because only single process is used now
Info setup(int NX, int NY, int NZ, int P, int STEPS) {
    //set openmp processors
    int processors_m = omp_get_num_procs();
    omp_set_num_threads(processors_m);
    //get paetition on every side
    int x_par;
    int y_par;
    int z_par;
    get_partition(NX, NY, NZ, P, &x_par, &y_par, &z_par, STEPS);
    Info result;
    result.x_par = x_par;
    result.y_par = y_par;
    result.z_par = z_par;
    int myrank = 0;
    // init
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int total_node = x_par * y_par * z_par;
    if (myrank < total_node){
        int node_x;
        int node_y;
        int node_z;
        int interval_x = (int) (NX / (float)x_par + 0.5) + 1;
        int interval_y = (int) (NY / (float)y_par + 0.5) + 1;
        int interval_z = (int) (NZ / (float)z_par + 0.5) + 1;
        int par_yz = y_par * z_par;
        node_x = myrank / (par_yz);
        int tmp = myrank % par_yz;
        node_y = tmp / (z_par);
        node_z = tmp % z_par;

        ll Asize_YZ = NY * NZ;

        result.offset_x = interval_x * node_x;
        result.end_x = min(NX, result.offset_x + interval_x);
        result.offset_y = interval_y * node_y;
        result.end_y = min(NY, result.offset_y + interval_y);
        result.offset_z = interval_z * node_z;
        result.end_z = min(NZ, result.offset_z + interval_z);
        result.nx = result.end_x - result.offset_x;
        result.ny = result.end_y - result.offset_y;
        result.nz = result.end_z - result.offset_z;


        //init the data area
        ll cube_size = (ll)result.nx * (ll)result.ny * (ll)result.nz;
        cube_blockA = (double *) malloc(cube_size * (ll)sizeof(double));
        cube_blockB = (double *) malloc(cube_size * (ll)sizeof(double));
        printf("cube size: %d \n", cube_size);
        int size_xy = result.nx * result.ny;
        int size_xz = result.nx * result.nz;
        int size_yz = result.ny * result.nz;
        int buffer_size = 2 * size_xy + 2 * size_xz + 2 * size_yz;
        receive_buffer = (double *) malloc(buffer_size * (ll) sizeof(double));
        printf("buffer size: %d \n", buffer_size);
        //send_buffer = (double *) malloc(buffer_size * (ll) sizeof(double));
        //begin to check the node to communicate
        if (result.offset_x != 0){
            node_rank[0] = myrank - par_yz;
        } else {
            node_rank[0] = -1;
        }
        if (result.end_x != NX){
            node_rank[2] = myrank + par_yz;
        } else {
            node_rank[2] = -1;
        }
        if (result.offset_y != 0){
            node_rank[3] = myrank - z_par;
        } else {
            node_rank[3] = -1;
        }
        if (result.end_y != NY){
            node_rank[1] = myrank + z_par;
        } else {
            node_rank[1] = -1;
        }
        if (result.offset_z != 0){
            node_rank[4] = myrank - 1;
        } else {
            node_rank[4] = -1;
        }
        if (result.end_z != NZ){
            node_rank[5] = myrank + 1;
        } else {
            node_rank[5] = -1;
        }
        buffer_offset[0] = 0;
        buffer_offset[1] = buffer_offset[0] + size_yz;
        buffer_offset[2] = buffer_offset[1] + size_xz;
        buffer_offset[3] = buffer_offset[2] + size_yz;
        buffer_offset[4] = buffer_offset[3] + size_xz;
        buffer_offset[5] = buffer_offset[4] + size_xy;
        buffer_offset[6] = buffer_offset[5] + size_xy;
    } else {
        result.end_x = 0;
        result.end_y = 0;
        result.end_z = 0;
        result.offset_x = 0;
        result.offset_y = 0;
        result.offset_z = 0;
        result.nx = 0;
        result.ny = 0;
        result.nz = 0;
    }
    return result;
}

//Must be re-written, including all the parameters
//compute the core mid part of the whole matrix

int stencil(Info info, int steps, int NX, int NY, int NZ) {
    static double x_buffer[3][1200 * 700];
    printf("begin \n");
    int hash_map[6] = {2, 3, 0, 1, 5, 4};
    int sub_s;

    int myrank;
    int size_xy = info.nx * info.ny;
    int size_xz = info.nx * info.nz;
    int size_yz = info.ny * info.nz;
    printf("cube info, xy : %d; xz : %d; yz : %d \n", size_xy, size_xz, size_yz);
    if (info.offset_x != 0){
        node_dsize[0] = size_yz;
    } else {
        node_dsize[0] = 0;
    }
    if (info.end_x != NX){
        node_dsize[2] = size_yz;
    } else {
        node_dsize[2] = 0;
    }
    if (info.offset_y != 0){
        node_dsize[3] = size_xz;
    } else {
        node_dsize[3] = 0;
    }
    if (info.end_y != NY){
        node_dsize[1] = size_xz;
    } else {
        node_dsize[1] = 0;
    }
    if (info.offset_z != 0){
        node_dsize[4] = size_xy;
    } else {
        node_dsize[4] = 0;
    }
    if (info.end_z != NZ){
        node_dsize[5] = size_xy;
    } else {
        node_dsize[5] = 0;
    }

    for (sub_s = 0; sub_s < 6; sub_s++){
        if (node_rank[sub_s] >= 0){
            printf("node rank %d, data size %d \n", node_rank[sub_s], node_dsize[sub_s]);
        }
    }

    int buffer_size = 2 * size_xy + 2 * size_xz + 2 * size_yz;
    double *send_buffer = (double *)malloc(buffer_size * sizeof(double));
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int i, j, k, s;
    //int size_yz = info.ny * info.nz;
    int nx = info.nx;
    int ny = info.ny;
    int nz = info.nz;
    for (s = 0; s < steps; s++) {
        //write data to send buffer and send
        //side 0
        //printf("begin to send data! \n");

        printf("-------------side 0-------------------\n");
        if (node_rank[0] >= 0){
            //printf("cal side 0! \n");
            double *start_index = (send_buffer + buffer_offset[0]);
            i = 0;
            int tmp_cube_i = i * size_yz;
            int tmp_cube_j;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_j, cube_index, count) schedule (dynamic)
            for (j = 0; j < ny; j++){
                tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_j + k;
                    count = cube_index - tmp_cube_i;
                    start_index[count] = cube_blockA[cube_index];
                }
            }
            MPI_Isend(start_index, node_dsize[0], MPI_DOUBLE, node_rank[0], (hash_map[0] + 1) * 10, MPI_COMM_WORLD, &send_request[0]);

        }
        // side 1
        printf("-------------side 1-------------------\n");
        if (node_rank[1] >= 0){
            //printf("cal side 1! \n");

            double *start_index = (send_buffer + buffer_offset[1]);
            j = ny - 1;
            int tmp_cube_j = j * nz;
            int tmp_cube_i;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_i, cube_index, count) schedule (dynamic)
            for (i = 0; i < nx; i++){
                tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_i + k;
                    count = i * nz + k;
                    start_index[count] = cube_blockA[cube_index];
                }
            }
            MPI_Isend(start_index, node_dsize[1], MPI_DOUBLE, node_rank[1], (hash_map[1] + 1) * 10, MPI_COMM_WORLD, &send_request[1]);

        }

        //side 2
        printf("-------------side 2-------------------\n");
        if (node_rank[2] >= 0){
            //printf("cal side 2! \n");

            double *start_index = (send_buffer + buffer_offset[2]);
            i = nx - 1;
            int tmp_cube_i = i * size_yz;
            int tmp_cube_j;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_j, cube_index, count) schedule (dynamic)
            for (j = 0; j < ny; j++){
                tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_j + k;
                    count = cube_index - tmp_cube_i;
                    start_index[count] = cube_blockA[cube_index];
                }
            }
            MPI_Isend(start_index, node_dsize[2], MPI_DOUBLE, node_rank[2], (hash_map[2] + 1) * 10, MPI_COMM_WORLD, &send_request[2]);

        }

        //side 3
        printf("-------------side 3-------------------\n");
        if (node_rank[3] >= 0){
            //printf("cal side 3! \n");
            double *start_index = (send_buffer + buffer_offset[3]);
            j = 0;
            int tmp_cube_j = j * nz;
            int tmp_cube_i;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_i, cube_index, count) schedule (dynamic)
            for (i = 0; i < nx; i++){
                tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_i + k;
                    count = i * nz + k;
                    start_index[count] = cube_blockA[cube_index];
                }
            }
            MPI_Isend(start_index, node_dsize[3], MPI_DOUBLE, node_rank[3], (hash_map[3] + 1) * 10, MPI_COMM_WORLD, &send_request[3]);
        }
        printf("-------------side 4-------------------\n");
        //side 4
        if (node_rank[4] >= 0){

            double *start_index = (send_buffer + buffer_offset[4]);
            k = 0;
            int tmp_cube_k = k;
            int tmp_cube_i;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_i, cube_index, count) schedule (dynamic)
            for (i = 0; i < nx; i++){
                tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++){
                    cube_index = tmp_cube_i + j * nz;
                    count = i * ny + j;
                    start_index[count] = cube_blockA[cube_index];
                }
                //printf("cal side 4, i : %d \n", i);
            }
            MPI_Isend(start_index, node_dsize[4], MPI_DOUBLE, node_rank[4], (hash_map[4] + 1) * 10, MPI_COMM_WORLD, &send_request[4]);
        }
        printf("-------------side 5-------------------\n");
        //side 5
        if (node_rank[5] >= 0) {

            double *start_index = (double *)(send_buffer + buffer_offset[5]);
            k = nz - 1;
            int tmp_cube_k = k;
            int tmp_cube_i;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_i, cube_index, count) schedule (dynamic)
            for (i = 0; i < nx; i++) {
                tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++) {
                    cube_index = tmp_cube_i + j * nz;
                    count = i * ny + j;
                    start_index[count] = cube_blockA[cube_index];
                }
            }
            MPI_Isend(start_index, node_dsize[5], MPI_DOUBLE, node_rank[5], (hash_map[5] + 1) * 10, MPI_COMM_WORLD, &send_request[5]);
        }
        // gather side data asyn
        int sub_s;
        for (sub_s = 0; sub_s < 6; sub_s++){
            if (node_rank[sub_s] >= 0){
                MPI_Irecv((receive_buffer + buffer_offset[sub_s]), node_dsize[sub_s], MPI_DOUBLE, node_rank[sub_s], (sub_s + 1) * 10, MPI_COMM_WORLD, &recv_request[sub_s]);
            }
        }


        printf("begin to compute \n");
        int total = omp_get_num_procs();

        #pragma omp parallel private(x_buffer)
        {

            int id = omp_get_thread_num();
            //printf("thread id : %d \n", id);
            int xh_interval = nx / (total / 2);
            int yz_no = id % 2;
            int zh_start;
            int zh_end;
            int yh_start;
            int yh_end;
            int xh_start = (id / 2) * xh_interval;
            int xh_end = min(xh_start + xh_interval, nx);
            if (nz < ny){
                zh_start = 0;
                zh_end = nz;
                if (yz_no == 0){
                    yh_start = 0;
                    yh_end = ny / 2 + 1;
                } else {
                    yh_start = ny / 2 - 1;
                    yh_end = ny;
                }
            } else {
                yh_start = 0;
                yh_end = ny;
                if (yz_no == 0){
                    zh_start = 0;
                    zh_end = nz / 2 + 1;
                } else {
                    zh_start = nz / 2 - 1;
                    zh_end = nz;
                }
            }
            //printf("thread id : %d, init done! \n", id);

            int pos_i;
            int count;
            for (pos_i = xh_start; pos_i < xh_end; pos_i++){
                int tmp_i = pos_i * size_yz;

                // begin to init cache
                if (pos_i != nx){
                    //plus 1
                    int x_buffer_index = (pos_i + 1) % 3;
                    int tmp_b_i = tmp_i + size_yz;
                    count = 0;
                    int pos_j;

                    for (pos_j = yh_start; pos_j < yh_end; pos_j++){
                        int tmp_b_j = tmp_b_i + pos_j * nz;
                        int pos_k;
                        for (pos_k = zh_start; pos_k < zh_end; pos_k++){
                            int tmp_b_index = tmp_b_j + pos_k;
                            //double tmp_sc = x_buffer[x_buffer_index][count];
                            //double a = tmp_sc + 1;
                            //printf("thread id : %d, x_buffer_index : %d, count : %d !\n", id, x_buffer_index, count);
                            x_buffer[x_buffer_index][count] = cube_blockA[tmp_b_index];
                            count++;
                        }
                    }
                }
                //printf("thread id : %d, init done_1! \n", id);

                if (pos_i == 0){
                    int tmp_b_i = 0;
                    int count = 0;
                    int pos_j;
                    for (pos_j = yh_start; pos_j < yh_end; pos_j++){
                        int tmp_b_j = tmp_b_i + pos_j * nz;
                        int pos_k;
                        for (pos_k = zh_start; pos_k < zh_end; pos_k++){
                            int tmp_b_index = tmp_b_j + pos_k;
                            x_buffer[0][count] = cube_blockA[tmp_b_index];
                            count++;
                        }
                    }
                }
                //printf("thread id : %d, init done_2! \n", id);
                //begin to cal
                int cache_index;
                int loop_y_start;
                int loop_y_end;
                if (yh_start == 0){
                    loop_y_start = yh_start;
                } else {
                    loop_y_start = yh_start + 1;
                }
                if (yh_end == ny){
                    loop_y_end = yh_end;
                } else {
                    loop_y_end = yh_end - 1;
                }
                int loop_z_start;
                if (zh_start == 0){
                    loop_z_start = zh_start;
                } else {
                    loop_z_start = zh_start + 1;
                }
                int loop_z_end;
                if (zh_end == nz){
                    loop_z_end = zh_end;
                } else {
                    loop_z_end = zh_end - 1;
                }
                int th_len_z = zh_end - zh_start;
                double sum;
                double r;
                int buffer_yz_index;
                int current_layer = pos_i % 3;
                int below_layer = (pos_i - 1) % 3;
                int above_layer = (pos_i + 1) % 3;
                int pos_j;
                //printf("thread id : %d, init done_3! \n", id);
                for (pos_j = loop_y_start; pos_j < loop_y_end; pos_j++){
                    ll tmp_j = tmp_i + pos_j * nz;
                    int pos_k;
                    for (pos_k = loop_z_start; pos_k < loop_z_end; pos_k++){
                        sum = 0;
                        ll tmp_index = tmp_j + pos_k;
                        buffer_yz_index = (pos_j - yh_start) * th_len_z + pos_k - zh_start;
                        r = x_buffer[current_layer][buffer_yz_index];
                        if (pos_k != nz - 1)
                            sum += x_buffer[current_layer][buffer_yz_index + 1];
                        if (pos_k != 0)
                            sum += x_buffer[current_layer][buffer_yz_index - 1];
                        // j + - 1
                        if (pos_j != ny - 1)
                            sum += x_buffer[current_layer][buffer_yz_index + th_len_z];
                        if (pos_j != 0)
                            sum += x_buffer[current_layer][buffer_yz_index - th_len_z];
                        // i + - 1
                        if (pos_i != nx - 1)
                            sum += x_buffer[above_layer][buffer_yz_index];
                        if (pos_i != 0)
                            sum += x_buffer[below_layer][buffer_yz_index];
                        cube_blockB[tmp_index] = r * 0.4 + sum * 0.1;
                    }
                }
                //printf("thread id : %d, init done_4! \n", id);
            }


        }

        /*
        ll tmp_i;
        ll tmp_j;
        ll tmp_index;
        #pragma omp parallel for private(tmp_i, tmp_j, tmp_index, r, sum) schedule (dynamic)
        for (i = 0; i < nx; i++){
            tmp_i = i * size_yz;
            for (j = 0; j < ny; j++){
                tmp_j = tmp_i + j * nz;
                for (k = 0; k < nz; k++){
                    tmp_index = tmp_j + k;
                    sum = 0;
                    r = cube_blockA[tmp_index];
                    // k + - 1
                    if (k != nz - 1)
                        sum += cube_blockA[tmp_index + 1];
                    if (k != 0)
                        sum += cube_blockA[tmp_index - 1];
                    // j + - 1
                    if (j != ny - 1)
                        sum += cube_blockA[tmp_index + nz];
                    if (j != 0)
                        sum += cube_blockA[tmp_index - nz];
                    // i + - 1
                    if (i != nx - 1)
                        sum += cube_blockA[tmp_index + size_yz];
                    if (i != 0)
                        sum += cube_blockA[tmp_index - size_yz];
                    cube_blockB[tmp_index] = r * 0.4 + sum * 0.1;
                }
            }
        }
         */
        //waite data to be received


        MPI_Status status;
        for (sub_s = 0; sub_s < 6; sub_s++){
            if (node_rank[sub_s] >= 0){
                //printf("waite receive %d \n", node_rank[sub_s]);
                MPI_Wait(&recv_request[sub_s], &status);

            }
        }
        // begin to cal the edge
        //side 0

        if (node_rank[0] >= 0){
            //printf("phase 2 cal side 0\n");

            double *start_index = (receive_buffer + buffer_offset[0]);
            i = 0;
            int tmp_cube_i = i * size_yz;
            int tmp_cube_j;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_j, cube_index, count) schedule (dynamic)
            for (j = 0; j < ny; j++){
                tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_j + k;
                    count = cube_index - tmp_cube_i;
                    cube_blockB[cube_index] += 0.1 * start_index[count];
                }
            }
        } else {
            i = 0;
            int tmp_cube_i = i * size_yz;
            int tmp_cube_j;
            int cube_index;
            #pragma omp parallel for private(tmp_cube_j, cube_index) schedule (dynamic)
            for (j = 0; j < ny; j++){
                tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_j + k;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }

        //side 1
        if (node_rank[1] >= 0){

            //printf("phase 2 cal side 1\n");

            double *start_index = (receive_buffer + buffer_offset[1]);
            //printf("receive_buffer[buffer_set[1]] : %lf \n", start_index[0]);
            j = ny - 1;
            int tmp_cube_j = j * nz;
            int tmp_cube_i;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_i, cube_index, count) schedule (dynamic)
            for (i = 0; i < nx; i++){
                tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_i + k;
                    count = i * nz + k;
                    cube_blockB[cube_index] += 0.1 * start_index[count];
                }
            }

        } else {
            j = ny - 1;
            int tmp_cube_j = j * nz;
            int tmp_cube_i;
            int cube_index;
            #pragma omp parallel for private(tmp_cube_i, cube_index) schedule (dynamic)
            for (i = 0; i < nx; i++){
                tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_i + k;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }

        //side 2
        if (node_rank[2] >= 0){
            //printf("phase 2 cal side 2\n");

            double *start_index = (receive_buffer + buffer_offset[2]);
            i = nx - 1;
            int tmp_cube_i = i * size_yz;
            int tmp_cube_j;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_j, cube_index, count) schedule (dynamic)
            for (j = 0; j < ny; j++){
                tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_j + k;
                    count = cube_index - tmp_cube_i;
                    cube_blockB[cube_index] += 0.1 * start_index[count];
                }
            }
        } else {
            i = nx - 1;
            int tmp_cube_i = i * size_yz;
            int tmp_cube_j;
            int cube_index;
            #pragma omp parallel for private(tmp_cube_j, cube_index) schedule (dynamic)
            for (j = 0; j < ny; j++){
                tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_j + k;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }
        //side 3
        if (node_rank[3] >= 0){

            double *start_index = (receive_buffer + buffer_offset[3]);
            j = 0;
            int tmp_cube_j = j * nz;
            int tmp_cube_i;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_i, cube_index, count) schedule (dynamic)
            for (i = 0; i < nx; i++){
                tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_i + k;
                    count = i * nz + k;
                    cube_blockB[cube_index] += 0.1 * start_index[count];
                }
            }
        } else {
            j = 0;
            int tmp_cube_j = j * nz;
            int tmp_cube_i;
            int cube_index;
            #pragma omp parallel for private(tmp_cube_j, cube_index) schedule (dynamic)
            for (i = 0; i < nx; i++){
                tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    cube_index = tmp_cube_i + k;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }

        //side 4
        if (node_rank[4] >= 0){
            double *start_index = (receive_buffer + buffer_offset[4]);
            k = 0;
            int tmp_cube_k = k;
            int tmp_cube_i;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_i, cube_index, count) schedule (dynamic)
            for (i = 0; i < nx; i++){
                tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++){
                    cube_index = tmp_cube_i + j * nz;
                    count = i * ny + j;
                    cube_blockB[cube_index] += 0.1 * start_index[count];
                }
            }
        } else {
            k = 0;
            int tmp_cube_k = k;
            int tmp_cube_i;
            int cube_index;
            #pragma omp parallel for private(tmp_cube_i, cube_index) schedule (dynamic)
            for (i = 0; i < nx; i++){
                tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++){
                    cube_index = tmp_cube_i + j * nz;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }
        //side 5
        if (node_rank[5] >= 0){
            //printf("phase 2 cal side 5\n");

            double *start_index = (receive_buffer + buffer_offset[5]);
            k = nz - 1;
            int tmp_cube_k = k;
            int tmp_cube_i;
            int cube_index;
            int count;
            #pragma omp parallel for private(tmp_cube_i, cube_index, count) schedule (dynamic)
            for (i = 0; i < nx; i++){
                tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++){
                    cube_index = tmp_cube_i + j * nz;
                    count = i * ny + j;
                    cube_blockB[cube_index] += 0.1 * start_index[count];
                }
            }
        } else {
            k = nz - 1;
            int tmp_cube_k = k;
            int tmp_cube_i;
            int cube_index;
            #pragma omp parallel for private(tmp_cube_i, cube_index) schedule (dynamic)
            for (i = 0; i < nx; i++){
                tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++){
                    cube_index = tmp_cube_i + j * nz;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }

        //check last send complete
        for (sub_s = 0; sub_s < 6; sub_s++){
            if (node_rank[sub_s] >= 0){
                MPI_Wait(&send_request[sub_s], &status);
            }
        }

        double *tmp = NULL;
        tmp = cube_blockA, cube_blockA = cube_blockB, cube_blockB = tmp;

    }
    free(send_buffer);
    return 0;
}

int main(int argc, char **argv) {
    double *send_buffer;// store the data to be send, 6 segments, not all are used
    int myrank, nprocs;
    int NX = 100, NY = 100, NZ = 100, STEPS = 10;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    //printf("mpi size %d \n", nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    NX = atoi(argv[1]), NY = atoi(argv[2]), NZ = atoi(argv[3]);
    STEPS = atoi(argv[4]);
    Info info = setup(NX, NY, NZ, nprocs, STEPS); //This is a single version
    if (myrank == 0){
        printf("partition info, x_par : %d; y_par : %d; z_par : %d \n", info.x_par, info.y_par, info.z_par);
    }
    total_nodes = info.x_par * info.y_par * info.z_par;
    //long long size = NX * NY * NZ;
    long long  size = info.nx * info.ny * info.nz;
    if (myrank == 0)
        printf("Size:%dx%dx%d, # of Steps: %d, # of procs: %d\n", NX, NY, NZ, STEPS, nprocs);
    cube_blockA = (double *) malloc(size * sizeof(double));
    Init(cube_blockA, size);
    printf("init data done! total nodes %d\n", total_nodes);
    struct timeval t1, t2;
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t1, NULL);
    printf("done 2 \n");
    if (myrank < total_nodes){
        printf("done 3 \n");
        stencil(info, STEPS, NX, NY, NZ);
        printf("done 4 \n");
    }
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t2, NULL);
    printf("barrier done! \n");
    if (myrank == 0) {
        printf("Total time: %.6lf\n", TIME(t1, t2));
    }
    Check(cube_blockA, size);
    free(cube_blockA);

    if (myrank < total_nodes){
        free(cube_blockB);
        free(receive_buffer);
    }

    printf("all done!\n");
    MPI_Finalize();
}
