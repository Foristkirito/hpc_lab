#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "./engine/proj_engine.h"


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

void copy_data(double *A, Info info, int NX, int NY, int NZ, int myrank){

    ll size_yz_T = NY * NZ;
    int size_yz = info.ny * info.nz;
    int i, j, k;
    int sub_i, sub_j, sub_k;
    printf("begin to copy on node %d \n", myrank);
    int count = 0;
    for (sub_i = info.offset_x; sub_i < info.end_x; sub_i++){
        //int tmp_cube_i = (sub_i - info.offset_x) * size_yz;
        ll tmp_A_i = sub_i * size_yz_T;
        for (sub_j = info.offset_y; sub_j < info.end_y; sub_j++){
            //int tmp_cube_j = tmp_cube_i + (sub_j - info.offset_y) * info.nz;
            ll tmp_A_j = tmp_A_i + sub_j * NZ;
            for (sub_k = info.offset_z; sub_k < info.end_z; sub_k++){
                //int cube_index = tmp_cube_j + sub_k - info.offset_z;
                ll A_index = tmp_A_j + sub_k;
                cube_blockA[count] = A[A_index];
                count++;
            }
        }
    }
    // copy data to six sides

    // side 0
    if (node_rank[0] >= 0){
        printf("copy side 0\n");
        int count = 0;
        double *start_index = (receive_buffer + buffer_offset[0]);
        sub_i = info.offset_x - 1;
        ll tmp_A_i = sub_i * size_yz_T;
        for (sub_j = info.offset_y; sub_j < info.end_y; sub_j++){
            ll tmp_A_j = tmp_A_i + sub_j * NZ;
            for (sub_k = info.offset_z; sub_k < info.end_z; sub_k++){
                ll A_index = tmp_A_j + sub_k;
                start_index[count] = A[A_index];
                count++;
            }
        }
        printf("receive[%d] : %lf \n", TESTN, start_index[TESTN]);
    }

    //side 1
    if (node_rank[1] >= 0){
        printf("copy side 1\n");
        int count = 0;
        double *start_index = (receive_buffer + buffer_offset[1]);
        //printf("done 1\n");
        sub_j = info.end_y;
        ll tmp_A_j = sub_j * NZ;
        //printf("offset x: %d, end_x: %d \n", info.offset_x, info.end_x);
        //printf("offset z: %d, end_z: %d \n", info.offset_z, info.end_z);
        for (sub_i = info.offset_x; sub_i < info.end_x; sub_i++){
            //printf("done 2_1\n");
            ll tmp_A_i = tmp_A_j + sub_i * size_yz_T;
            for (sub_k = info.offset_z; sub_k < info.end_z; sub_k++){
                ll A_index = tmp_A_i + sub_k;
                //printf("done 2\n");
                start_index[count] = A[A_index];
                //printf("A[A_index] : %lf; start_index[%d] : %lf \n", A[A_index], count, start_index[count]);
                count++;
            }
        }
        printf("receive[%d] : %lf \n", TESTN, start_index[TESTN]);
    }

    //side 2
    if (node_rank[2] >= 0){
        printf("copy side 2\n");
        int count = 0;
        double *start_index = (receive_buffer + buffer_offset[2]);
        sub_i = info.end_x;
        ll tmp_A_i = sub_i * size_yz_T;
        for (sub_j = info.offset_y; sub_j < info.end_y; sub_j++){
            ll tmp_A_j = tmp_A_i + sub_j * NZ;
            for (sub_k = info.offset_z; sub_k < info.end_z; sub_k++){
                ll A_index = tmp_A_j + sub_k;
                start_index[count] = A[A_index];
                count++;
            }
        }
        printf("receive[%d] : %lf \n", TESTN, start_index[TESTN]);
    }

    //side 3
    if (node_rank[3] >= 0){
        printf("copy side 3\n");
        int count = 0;
        double *start_index = (receive_buffer + buffer_offset[3]);
        sub_j = info.offset_y - 1;
        ll tmp_A_j = sub_j * NZ;
        for (sub_i = info.offset_x; sub_i < info.end_x; sub_i++){
            ll tmp_A_i = tmp_A_j + sub_i * size_yz_T;
            for (sub_k = info.offset_z; sub_k < info.end_z; sub_k++){
                ll A_index = tmp_A_i + sub_k;
                start_index[count] = A[A_index];
                //printf("A[A_index] : %lf; start_index[%d] : %lf \n", A[A_index], count, start_index[count]);
                count++;
            }
        }
        printf("receive[%d] : %lf \n", TESTN, start_index[TESTN]);
    }

    //side 4
    if (node_rank[4] >= 0){
        printf("copy side 4\n");
        int count = 0;
        double *start_index = (receive_buffer + buffer_offset[4]);
        sub_k = info.offset_z - 1;
        ll tmp_A_k = sub_k;
        for (sub_i = info.offset_x; sub_i < info.end_x; sub_i++){
            ll tmp_A_i = tmp_A_k + sub_i * size_yz_T;
            for (sub_j = info.offset_y; sub_j < info.end_y; sub_j++){
                ll A_index = tmp_A_i + sub_j * NZ;
                start_index[count] = A[A_index];
                count++;
            }
        }
    }

    // side 5
    if (node_rank[5] >= 0){
        printf("copy side 5\n");
        int count = 0;
        double *start_index = (receive_buffer + buffer_offset[5]);
        sub_k = info.end_z;
        ll tmp_A_k = sub_k;
        for (sub_i = info.offset_x; sub_i < info.end_x; sub_i++){
            ll tmp_A_i = tmp_A_k + sub_i * size_yz_T;
            for (sub_j = info.offset_y; sub_j < info.end_y; sub_j++){
                ll A_index = tmp_A_i + sub_j * NZ;
                start_index[count] = A[A_index];
                count++;
            }
        }
    }
}

//Must be re-written, including all the parameters
//compute the core mid part of the whole matrix

int stencil(double *A, Info info, int steps, int NX, int NY, int NZ) {
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
        printf("step %d \n", s);
        if (s != 0){
            //before cal get the data asyn

            int sub_s;
            for (sub_s = 0; sub_s < 6; sub_s++){
                if (node_rank[sub_s] >= 0){
                    //printf("receive size %d \n", node_dsize[sub_s]);
                    MPI_Irecv((receive_buffer + buffer_offset[sub_s]), node_dsize[sub_s], MPI_DOUBLE, node_rank[sub_s], (sub_s + 1) * 10, MPI_COMM_WORLD, &recv_request[sub_s]);
                }
            }
        }
        //printf("begin to cal phase_1! \n");
        for (i = 0; i < nx; i++){
            ll tmp_i = i * size_yz;
            for (j = 0; j < ny; j++){
                ll tmp_j = tmp_i + j * nz;
                for (k = 0; k < nz; k++){
                    ll tmp_index = tmp_j + k;
                    double r = 0.4 * cube_blockA[tmp_index];
                    // k + - 1
                    if (k != nz - 1)
                        r += 0.1 * cube_blockA[tmp_index + 1];
                    if (k != 0)
                        r += 0.1 * cube_blockA[tmp_index - 1];
                    // j + - 1
                    if (j != ny - 1)
                        r += 0.1 * cube_blockA[tmp_index + nz];
                    if (j != 0)
                        r += 0.1 * cube_blockA[tmp_index - nz];
                    // i + - 1
                    if (i != nx - 1)
                        r += 0.1 * cube_blockA[tmp_index + size_yz];
                    if (i != 0)
                        r += 0.1 * cube_blockA[tmp_index - size_yz];
                    cube_blockB[tmp_index] = r;
                }
            }
        }
        if (s != 0){
            //printf("begin to wait receive! \n");
            //sleep(60000);
            //ensure all the edge needed has been received
            int sub_s;
            MPI_Status status;
            for (sub_s = 0; sub_s < 6; sub_s++){
                if (node_rank[sub_s] >= 0){
                    //printf("waite receive %d \n", node_rank[sub_s]);
                    MPI_Wait(&recv_request[sub_s], &status);

                }
            }
        }
       // printf("begin to cal phase_2! \n");
        // begin to cal the edge
        //side 0
        if (node_rank[0] >= 0){
            //printf("phase 2 cal side 0\n");
            int count = 0;
            double *start_index = (receive_buffer + buffer_offset[0]);
            i = 0;
            int tmp_cube_i = i * size_yz;
            for (j = 0; j < ny; j++){
                int tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_j + k;
                    double pre = cube_blockB[cube_index];
                    cube_blockB[cube_index] += 0.1 * start_index[count];
                    count++;
                }
            }
        } else {
            i = 0;
            int tmp_cube_i = i * size_yz;
            for (j = 0; j < ny; j++){
                int tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_j + k;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }

        //side 1
        if (node_rank[1] >= 0){

            //printf("phase 2 cal side 1\n");
            int count = 0;
            double *start_index = (receive_buffer + buffer_offset[1]);
            //printf("receive_buffer[buffer_set[1]] : %lf \n", start_index[0]);
            j = ny - 1;
            int tmp_cube_j = j * nz;
            for (i = 0; i < nx; i++){
                int tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_i + k;
                    //printf("start_index[count] : %lf\n", start_index[count]);
                    double pre = cube_blockB[cube_index];
                    cube_blockB[cube_index] = pre + 0.1 * start_index[count];
                    count++;
                }
            }

        } else {
            j = ny - 1;
            int tmp_cube_j = j * nz;
            for (i = 0; i < nx; i++){
                int tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_i + k;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }

        //side 2
        if (node_rank[2] >= 0){
            //printf("phase 2 cal side 2\n");
            int count = 0;
            double *start_index = (receive_buffer + buffer_offset[2]);
            i = nx - 1;
            int tmp_cube_i = i * size_yz;
            for (j = 0; j < ny; j++){
                int tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_j + k;
                    double pre = cube_blockB[cube_index];
                    cube_blockB[cube_index] += 0.1 * start_index[count];
                    count++;
                }
            }
        } else {
            i = nx - 1;
            int tmp_cube_i = i * size_yz;
            for (j = 0; j < ny; j++){
                int tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_j + k;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }
        //side 3
        if (node_rank[3] >= 0){

            //printf("phase 2 cal side 3\n");
            int count = 0;
            double *start_index = (receive_buffer + buffer_offset[3]);
            j = 0;
            int tmp_cube_j = j * nz;
            for (i = 0; i < nx; i++){
                int tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_i + k;
                    double pre = cube_blockB[cube_index];
                    cube_blockB[cube_index] = pre + 0.1 * start_index[count];
                    //printf("pre_cube : %lf; post_cube: %lf \n", pre, cube_blockB[cube_index]);
                    count++;
                }
            }


        } else {
            j = 0;
            int tmp_cube_j = j * nz;
            for (i = 0; i < nx; i++){
                int tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_i + k;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }

        //side 4
        if (node_rank[4] >= 0){
            //printf("phase 2 cal side 4\n");
            int count = 0;
            double *start_index = (receive_buffer + buffer_offset[4]);
            k = 0;
            int tmp_cube_k = k;
            for (i = 0; i < nx; i++){
                int tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++){
                    int cube_index = tmp_cube_i + j * nz;
                    cube_blockB[cube_index] += 0.1 * start_index[count];
                    count++;
                }
            }
        } else {
            k = 0;
            int tmp_cube_k = k;
            for (i = 0; i < nx; i++){
                int tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++){
                    int cube_index = tmp_cube_i + j * nz;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }
        //side 5
        if (node_rank[5] >= 0){
            //printf("phase 2 cal side 5\n");
            int count = 0;
            double *start_index = (receive_buffer + buffer_offset[5]);
            k = nz - 1;
            int tmp_cube_k = k;
            for (i = 0; i < nx; i++){
                int tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++){
                    int cube_index = tmp_cube_i + j * nz;
                    cube_blockB[cube_index] += 0.1 * start_index[count];
                    count++;
                }
            }
        } else {
            k = nz - 1;
            int tmp_cube_k = k;
            for (i = 0; i < nx; i++){
                int tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++){
                    int cube_index = tmp_cube_i + j * nz;
                    cube_blockB[cube_index] += 0.1 * cube_blockA[cube_index];
                }
            }
        }
        //check last send complete
        if (s != 0){
            //printf("begin to wait last send done! \n");
            int sub_s;
            MPI_Status status;
            for (sub_s = 0; sub_s < 6; sub_s++){
                if (node_rank[sub_s] >= 0){
                    MPI_Wait(&send_request[sub_s], &status);
                }
            }
        }
        //write data to send buffer and send
        //side 0
        //printf("begin to send data! \n");
        //printf("-------------side 0-------------------\n");
        if (node_rank[0] >= 0){
            //printf("cal side 0! \n");
            int count = 0;
            double *start_index = (send_buffer + buffer_offset[0]);
            i = 0;
            int tmp_cube_i = i * size_yz;
            for (j = 0; j < ny; j++){
                int tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_j + k;
                    start_index[count] = cube_blockB[cube_index];
                    count++;
                }
            }
            if (s != steps - 1){
                //printf("send side 0! %d data to send \n", node_dsize[0]);
                //send
                MPI_Isend(start_index, node_dsize[0], MPI_DOUBLE, node_rank[0], (hash_map[0] + 1) * 10, MPI_COMM_WORLD, &send_request[0]);
            }

        }
        // side 1
        //printf("-------------side 1-------------------\n");
        if (node_rank[1] >= 0){
            //printf("cal side 1! \n");
            int count = 0;
            double *start_index = (send_buffer + buffer_offset[1]);
            j = ny - 1;
            int tmp_cube_j = j * nz;
            for (i = 0; i < nx; i++){
                int tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_i + k;
                    start_index[count] = cube_blockB[cube_index];
                    count++;
                }
            }
            if (s != steps - 1){
                //send
                //printf("send side 1! \n");
                MPI_Isend(start_index, node_dsize[1], MPI_DOUBLE, node_rank[1], (hash_map[1] + 1) * 10, MPI_COMM_WORLD, &send_request[1]);
            }

        }

        //side 2
        //printf("-------------side 2-------------------\n");
        if (node_rank[2] >= 0){
            //printf("cal side 2! \n");
            int count = 0;
            double *start_index = (send_buffer + buffer_offset[2]);
            i = nx - 1;
            int tmp_cube_i = i * size_yz;
            for (j = 0; j < ny; j++){
                int tmp_cube_j = tmp_cube_i + j * nz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_j + k;
                    start_index[count] = cube_blockB[cube_index];
                    count++;
                }
            }
            if (s != steps - 1){
                //printf("send side 2! %d to send! \n", node_dsize[2]);
                //send
                MPI_Isend(start_index, node_dsize[2], MPI_DOUBLE, node_rank[2], (hash_map[2] + 1) * 10, MPI_COMM_WORLD, &send_request[2]);
            }

        }

        //side 3
        //printf("-------------side 3-------------------\n");
        if (node_rank[3] >= 0){
            //printf("cal side 3! \n");
            int count = 0;
            double *start_index = (send_buffer + buffer_offset[3]);
            j = 0;
            int tmp_cube_j = j * nz;
            for (i = 0; i < nx; i++){
                int tmp_cube_i = tmp_cube_j + i * size_yz;
                for (k = 0; k < nz; k++){
                    int cube_index = tmp_cube_i + k;
                    start_index[count] = cube_blockB[cube_index];
                    count++;
                }
            }
            if (s != steps - 1){
                //printf("send side 3! \n");
                //send
                MPI_Isend(start_index, node_dsize[3], MPI_DOUBLE, node_rank[3], (hash_map[3] + 1) * 10, MPI_COMM_WORLD, &send_request[3]);
            }
        }
        //printf("-------------side 4-------------------\n");
        //side 4
        if (node_rank[4] >= 0){
            int count = 0;
            double *start_index = (send_buffer + buffer_offset[4]);
            k = 0;
            int tmp_cube_k = k;
            //printf("cal side 4! \n");
            for (i = 0; i < nx; i++){
                int tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++){
                    int cube_index = tmp_cube_i + j * nz;
                    start_index[count] = cube_blockB[cube_index];
                    count++;
                }
                //printf("cal side 4, i : %d \n", i);
            }
            if (s != steps - 1){
                //printf("send side 4! \n");
                //send
                MPI_Isend(start_index, node_dsize[4], MPI_DOUBLE, node_rank[4], (hash_map[4] + 1) * 10, MPI_COMM_WORLD, &send_request[4]);
            }
        }
        //printf("-------------side 5-------------------\n");
        //side 5
        if (node_rank[5] >= 0) {
            int count = 0;
            double *start_index = (double *)(send_buffer + buffer_offset[5]);
            k = nz - 1;
            int tmp_cube_k = k;
            for (i = 0; i < nx; i++) {
                int tmp_cube_i = tmp_cube_k + i * size_yz;
                for (j = 0; j < ny; j++) {
                    int cube_index = tmp_cube_i + j * nz;
                    start_index[count] = cube_blockB[cube_index];
                    count++;
                }
            }
            if (s != steps - 1){
                //printf("send side 5! \n");
                //send
                MPI_Isend(start_index, node_dsize[5], MPI_DOUBLE, node_rank[5], (hash_map[5] + 1) * 10, MPI_COMM_WORLD, &send_request[5]);
            }
        }
        //printf("One turn done! \n");
        double *tmp = NULL;
        tmp = cube_blockA, cube_blockA = cube_blockB, cube_blockB = tmp;

    }
    free(send_buffer);
    return 0;
}

//node 0 get data from other node
void gather_data(Info info, double *A, int NX, int NY, int NZ){
    int x_par = info.x_par;
    int y_par = info.y_par;
    int z_par = info.z_par;
    int interval_x = (int) (NX / (float)x_par + 0.5) + 1;
    int interval_y = (int) (NY / (float)y_par + 0.5) + 1;
    int interval_z = (int) (NZ / (float)z_par + 0.5) + 1;
    int buffer_size = interval_x * interval_y * interval_z;
    double *get_buffer = (double *)malloc(buffer_size * sizeof(double));
    ll size_yz_T = NY * NZ;
    int size_yz = info.ny * info.nz;
    int i, j, k;
    int sub_i, sub_j, sub_k;
    printf("rank 0 begin to gather data \n");
    for (i = 0; i < x_par; i++){
        for (j = 0; j < y_par; j++){
            for (k = 0; k < z_par; k++){
                int t_node = i * y_par * z_par + j * z_par + k;
                //printf("target node %d \n", t_node);
                if (t_node != 0){
                    //use mpi to receive data
                    MPI_Status m_status;
                    MPI_Recv(get_buffer, buffer_size, MPI_DOUBLE, t_node, t_node, MPI_COMM_WORLD, &m_status);
                    //printf("data geted!   begin to write buffer \n");
                    int offset_x = i * interval_x;
                    int end_x = min(offset_x + interval_x, NX);
                    int offset_y = j * interval_y;
                    int end_y = min(offset_y + interval_y, NY);
                    int offset_z = k * interval_z;
                    int end_z = min(offset_z + interval_z, NZ);
                    int buffer_size = (end_x - offset_x) * (end_y - offset_y) * (end_z - offset_z);
                    int tmp_yz = (end_y - offset_y) * (end_z - offset_z);
                    int tmp_z = end_z - offset_z;
                    int count = 0;
                    for (sub_i = offset_x; sub_i < end_x; sub_i++){
                        ll tmp_i = sub_i * size_yz_T;
                        for (sub_j = offset_y; sub_j < end_y; sub_j++){
                            ll tmp_j = tmp_i + sub_j * NZ;
                            for (sub_k = offset_z; sub_k < end_z; sub_k++){
                                ll A_index = tmp_j + sub_k;
                                A[A_index] = get_buffer[count];
                                count++;
                            }
                        }
                    }
                } else {
                    int count = 0;
                    for (sub_i = info.offset_x; sub_i < info.end_x; sub_i++){
                        ll tmp_i = sub_i * size_yz_T;
                        for (sub_j = info.offset_y; sub_j < info.end_y; sub_j++){
                            ll tmp_j = tmp_i + sub_j * NZ;
                            for (sub_k = info.offset_z; sub_k < info.end_z; sub_k++){
                                ll A_index = tmp_j + sub_k;
                                A[A_index] = cube_blockA[count];
                                count++;
                            }
                        }
                    }
                }
            }
        }
    }
    free(get_buffer);
}

//other nodes begin to send data to rank 0;
void send_data(Info info, int myrank){
    printf("begin to send data to rank 0 \n");
    int cube_size = info.nx * info.ny * info.nz;
    MPI_Send(cube_blockA, cube_size, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
    printf("date transfered done! \n");
}

int main(int argc, char **argv) {
    double *A = NULL;
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
    long long size = NX * NY * NZ;
    if (myrank == 0)
        printf("Size:%dx%dx%d, # of Steps: %d, # of procs: %d\n", NX, NY, NZ, STEPS, nprocs);
    A = (double *) malloc(size * sizeof(double));
    Init(A, size);
    printf("init data done! total nodes %d\n", total_nodes);
    if (myrank < total_nodes){
        copy_data(A, info, NX, NY, NZ, myrank);
        printf("copy data done! \n");
    }
    struct timeval t1, t2;
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t1, NULL);
    if (myrank < total_nodes){
        stencil(A, info, STEPS, NX, NY, NZ);
    }
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t2, NULL);
    printf("barrier done! \n");
    if (myrank == 0) {
        printf("Total time: %.6lf\n", TIME(t1, t2));
    }
    if (myrank == 0){
        gather_data(info, A, NX, NY, NZ);
        printf("rank 0 gather data done!\n");
    } else {
        if (myrank < total_nodes){
            send_data(info, myrank);
            printf("send data done!\n");
        }
    }
    Check(A, size);
    free(A);

    if (myrank < total_nodes){
        free(cube_blockA);
        free(cube_blockB);
        free(receive_buffer);
    }

    printf("all done!\n");
    MPI_Finalize();
}
