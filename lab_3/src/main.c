#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/fcntl.h>
#include <string.h>
#include <mpi.h>

#define TIME(a, b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))
#define RST 0.000067247278541769


const int VECTOR_SIZE = 4;
const int k = 5;
typedef double vec
    __attribute__ ((vector_size (sizeof(double) * VECTOR_SIZE)));
typedef long long ll;
typedef struct {
  int myrank;
  int nx;
  int nz;
  int ny;
  int tX;
  int tY;
  int tZ;
  int x_par;
  int y_par;
  int z_par;
  int x_start;
  int x_end;
  int y_start;
  int y_end;
  int z_start;
  int z_end;
  int len_vb;
  int len_x;
} Data_Info;


typedef struct{
  int col;
  double val;
} pair;



double *A_value; // partial A
int *A_col;
int *A_ptr;
double *b; // partial b
double *x; // partial x, recieve buffer is behind the local x
double *send_buffer;//store the data to send, 8 segments, not all are used
MPI_Request recv_request[8]; // receive request to check needed data has been received. not all are used
MPI_Request send_request[8]; // send requets
int buffer_offset[8]; // the start index of par i to send
int node_dsize[8]; // this node communicates with 8 nodes at most, this array save the date size to transfer
int node_rank[8]; // save the node rank to communication
int target_side[8]; // 初始化目标的index, 从接受者的角度来看
// receive buffer is stored in R_hat

int x_cor[19] = {0, -1, +1, 0, 0, +1, +1, -1, -1, 0, -1, +1, 0, 0, 0, -1, +1, 0, 0};
int y_cor[19] = {0, 0, 0, -1, +1, +1, -1, -1, +1, 0, 0, 0, -1, +1, 0, 0, 0, -1, +1};
int z_cor[19] = {0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, +1, +1, +1, +1, +1};

int min(int a, int b){
    if (a > b){
        return b;
    } else {
        return a;
    }
}

Data_Info init_info(int NX, int NY, int NZ, int PX, int PY, int PZ){
    Data_Info info;
    info.tX = NX;
    info.tY = NY;
    info.tZ = NZ;
    info.x_par = PX;
    info.y_par = PY;
    info.z_par = PZ;
    info.nz = NZ;
    info.z_start = 0;
    info.z_end = NZ;
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    info.myrank = myrank;
    int sub_i;
    if (PX == 1 && PY == 1 && PZ == 1){
        //info.myrank
        info.nx = NX;
        info.ny = NY;
        info.x_start = 0;
        info.x_end = NX;
        info.y_start = 0;
        info.y_end = NY;
        ll size_col = (ll) NX * (ll) NY * (ll) NZ;
        ll size_a = size_col * 19;
        A_value = (double *) malloc(size_a * sizeof(double));
        A_col = (int *) malloc(size_a * sizeof(int));
        A_ptr = (int *) malloc((size_col + 1) * sizeof(int));
        b = (double *) malloc(size_col * sizeof(double));
        info.len_vb = (int)size_col;
        // cal buffer size
        int size_buffer = 4 * NZ + 2 * NZ * NY + 2 * NZ * NX;
        //receive_buffer = (double *) malloc(size_buffer * sizeof(double));
        info.len_x = info.len_vb + size_buffer;
        g_map = (int *) malloc(size_buffer * sizeof(int));
        for (sub_i = 0; sub_i < size_buffer; sub_i++){
            g_map[sub_i] = -1;
        }
        x = (double *) malloc(info.len_x * sizeof(double));
        memset(x, 0, info.len_x * sizeof(double));
        send_buffer = (double *) malloc(size_buffer * sizeof(double));
        //init node info
        int i;
        for (i = 0; i < 8; i++){
            node_dsize[i] = -1;
            node_rank[i] = 0;
            buffer_offset[i] = -1;
        }
        return info;
    } else {
        // 首先切 Y 轴, naive
        int par_y_index = myrank / PX;
        int interval_y = NY / PY + 1;
        info.y_start = par_y_index * interval_y;
        info.y_end = min(NY, info.y_start + interval_y);
        info.ny = info.y_end - info.y_start;
        // 切 X 轴, 比较繁琐
        int par_x_index = myrank % PX;
        int bi_par_x = PX / 2;
        int bi_nx = NX / 2;
        int interval_x = bi_nx / bi_par_x;
        if (par_x_index < bi_par_x){
            info.x_start = interval_x * par_x_index;
            info.x_end = min(bi_nx, info.x_start + interval_x);
        } else {
            info.x_start = bi_nx + interval_x * (par_x_index - bi_par_x);
            info.x_end = min(NX, info.x_start + interval_x);
            //printf("x_start: %d, x_end : %d \n", info.x_start, info.x_end);
        }
        info.nx = info.x_end - info.x_start;
    }
    int size_col = info.nx * info.ny * info.nz;
    int size_a = size_col * 19;
    A_value = (double *) malloc(size_a * sizeof(double));
    A_col = (int *) malloc(size_a * sizeof(int));
    A_ptr = (int *) malloc((size_col + 1) * sizeof(int));
    //printf("nx : %d, ny : %d, nz : %d\n", info.nx, info.ny, info.nz);
    b = (double *) malloc(size_col * sizeof(double));
    info.len_vb = (int)size_col;
    // cal buffer size
    int size_buffer = 4 * info.nz + 2 * info.nz * info.ny + 2 * info.nz * info.nx;
    g_map = (int *) malloc(size_buffer * sizeof(int));
    for (sub_i = 0; sub_i < size_buffer; sub_i++){
        g_map[sub_i] = -1;
    }
    if (g_map == NULL){
        printf("null\n");
    }
    //receive_buffer = (double *) malloc(size_buffer * sizeof(double));
    info.len_x = info.len_vb + size_buffer;
    x = (double *) malloc(info.len_x * sizeof(double));
    for (sub_i = 0; sub_i < info.len_x; sub_i++){
        x[sub_i] = -2;
    }
    memset(x, 0, info.len_x * sizeof(double));
    send_buffer = (double *) malloc(size_buffer * sizeof(double));
    //init buffer start index and buffer size, same to send and receive
    buffer_offset[0] = 0;
    buffer_offset[1] = info.nz + buffer_offset[0];
    buffer_offset[2] = info.nz * info.ny + buffer_offset[1];
    buffer_offset[3] = info.nz + buffer_offset[2];
    buffer_offset[4] = info.nz * info.nx + buffer_offset[3];
    buffer_offset[5] = info.nz + buffer_offset[4];
    buffer_offset[6] = info.nz * info.ny + buffer_offset[5];
    buffer_offset[7] = info.nz + buffer_offset[6];
    int i;
    for(i = 0; i < 7; i++){
        node_dsize[i] = buffer_offset[i + 1] - buffer_offset[i];
        //printf("node data size[%d] : %d \n", i, node_dsize[i]);
    }
    node_dsize[7] = info.nz * info.nx;
    // init node to send and receive, and the target side from the receiver
    int line_rank = myrank % PX;
    int line_lay = myrank / PX;
    int bi_px = PX / 2;
    if (info.y_start == 0){
        //side 0
        node_rank[0] = (line_rank - 1 + bi_px) % PX + line_lay * PX;
        target_side[0] = 6; //send
        // side 7
        node_rank[7] = (line_rank + bi_px) % PX + line_lay * PX;
        target_side[7] = 7;
        //side 6
        node_rank[6] = (line_rank + 1 + bi_px) % PX + line_lay * PX;
        target_side[6] = 0;
    } else {
        // side 0
        node_rank[0] = (line_rank - 1 + PX) % PX + line_lay * PX - PX;
        target_side[0] = 4;
        // side 7
        node_rank[7] = myrank - PX;
        target_side[7] = 3;
        //side 6
        node_rank[6] = (line_rank + 1) % PX + line_lay * PX - PX;
        target_side[6] = 2;
    }
    // side 1
    node_rank[1] = (line_rank - 1 + PX) % PX + line_lay * PX;
    target_side[1] = 5;
    // side 5
    node_rank[5] = (line_rank + 1) % PX + line_lay * PX;
    target_side[5] = 1;
    if (info.y_end == NY){
        //side 2
        node_rank[2] = (line_rank - 1 + bi_px) % PX + line_lay * PX;
        target_side[2] = 4; //send
        // side 3
        node_rank[3] = (line_rank + bi_px ) % PX + line_lay * PX;
        target_side[3] = 3;
        // side 4
        node_rank[4] = (line_rank + 1 + bi_px) % PX + line_lay * PX;
        target_side[4] = 2;
    } else {
        //side 2
        node_rank[2] = (line_rank - 1 + PX) % PX + line_lay * PX + PX;
        target_side[2] = 6;
        //side 3
        node_rank[3] = myrank + PX;
        target_side[3] = 7;
        //side 4
        node_rank[4] = (line_rank + 1) % PX + line_lay * PX + PX;
        target_side[4] = 0;
    }
    return info;
}

int get_col(int i, int j, int k, Data_Info info){
    if (i >= 0 && i < info.nx && j >= 0 && j < info.ny && k >=0 && k < info.nz){
        return (i * info.ny * info.nz + j * info.nz + k);
    }
    int g_i = i + info.x_start;
    int g_j = j + info.y_start;
    int g_k = k + info.z_start;
    if (g_k < 0 || g_k >= info.tZ){
        return -1;
    }
    if (g_j == -1){
        g_j = 0;
        g_i += (info.tX >> 1);
    } else {
        if (g_j == info.tY){
            g_j = info.tY - 1;
            g_i += (info.tX >> 1);
        }
    }
    // x 循环回来
    g_i = (g_i + info.tX) % info.tX;
    if (g_i >= info.x_start && g_i < info.x_end && g_j >= info.y_start && g_j < info.y_end && g_k >= info.z_start && g_k < info.z_end){
        // the node is still in cube
        int l_i = g_i - info.x_start;
        int l_j = g_j - info.y_start;
        int l_k = g_k - info.z_start;
        return (l_i * (info.ny * info.nz) + l_j * info.nz + l_k);
    } else {
        // the node is outof cube
        int size = info.nx * info.ny * info.nz;
        if (j == -1){
            if (i == -1){
                //side 0
                return (size + k);
            } else {
                if (i == info.nx){
                    //side 6
                    return (size + info.nz * 3 + info.nz * info.ny * 2 + info.nz * info.nx + k);
                } else {
                    //side 7
                    return (size + info.nz * 4 + info.nz * info.ny * 2 + info.nz * info.nx + i * info.nz + k);
                }
            }
        } else {
            if (j == info.ny){
                if (i == -1){
                    //side 2
                    return (size + info.nz + info.nz * info.ny + k);
                } else {
                    if (i == info.nx){
                        //side 4
                        return (size + info.nz * 2 + info.nz * info.ny + info.nz * info.nx + k);
                    } else {
                        //side 3
                        return (size + info.nz * 2 + info.nz * info.ny + i * info.nz + k);
                    }
                }
            } else {
                if (i == -1){
                    //side 1
                    return (size + info.nz + j * info.nz + k);
                } else {
                    //side 5
                    if (i == info.nx){
                        return (size + info.nz * 3 + info.nz * info.ny + info.nz * info.nx + j * info.nz + k);
                    } else {
                        //printf("error!!!!!!!!!!!!!!!!!!!!\n");
                        return -1;
                    }

                }
            }
        }
    }
}

int pair_comparitor (const void* lhs, const void* rhs) {
    return (((pair *)lhs)->col - ((pair *)rhs)->col);
}

void init_A(char *filename_A, Data_Info info){
    int input_fd = open(filename_A, O_RDONLY);
    ll size = info.nx * info.ny * info.nz;
    double *data = (double *)malloc(size * 19 * sizeof(double));
    ll file_offset = info.x_start * info.tY * info.nz + info.y_start * info.nz + info.z_start;
    ll result = lseek(input_fd, file_offset * 19 * sizeof(double), SEEK_SET);
    printf("file offset : %d\n", result);
    A_ptr[0] = 0;
    if ((read(input_fd, data, sizeof(double) * size * 19)) > 0){
        int i, j, k;
        ll size_cube = info.ny * info.nz;
        int count = 0;
        for (i = 0; i < info.nx; i++){
            for (j = 0; j < info.ny; j++){
                for (k = 0; k < info.nz; k++){
                    pair data_tmp[19];
                    int pair_count = 0;
                    ll index_cube = i * size_cube + j * info.nz + k; //row number of cube
                    //get col
                    int sub_num;
                    for (sub_num = 0; sub_num < 19; sub_num++){
                        double val = data[index_cube * 19 + sub_num];
                        int col = get_col(i + x_cor[sub_num], j + y_cor[sub_num], k + z_cor[sub_num], info);

                        if (col >= 0){
                            data_tmp[pair_count].val = val;
                            data_tmp[pair_count].col = col;
                            pair_count++;
                        }
                    }
                    // begin to order the pair
                    qsort(data_tmp, pair_count, sizeof(pair), pair_comparitor);
                    //begin to fill into A

                    for (sub_num = 0; sub_num < pair_count; sub_num++){
                        A_value[count] = data_tmp[sub_num].val;
                        A_col[count] = data_tmp[sub_num].col;
                        count++;
                    }
                    A_ptr[index_cube + 1] = A_ptr[index_cube] + pair_count;
                }
            }
        }
    } else {
        printf("read A file error!\n");
    }
    free(data);
    close(input_fd);
}

void init_b(char *filename_b, Data_Info info){
    int input_fd = open(filename_b, O_RDONLY);
    ll size = info.nx * info.ny * info.nz;
    double *data = (double *)malloc(size * sizeof(double));
    ll file_offset = info.x_start * info.tY * info.nz + info.y_start * info.nz + info.z_start;
    ll result = lseek(input_fd, file_offset * sizeof(double), SEEK_SET);
    if ((read(input_fd, data, sizeof(double) * size)) > 0){
        int i, j, k;
        ll size_cube = info.ny * info.nz;
        for (i = 0; i < info.nx; i++){
            for (j = 0; j < info.ny; j++){
                for (k = 0; k < info.nz; k++){
                    ll index_cube = i * size_cube + j * info.nz + k; //row number of cube
                    //fill into b in x y z order
                    b[index_cube] = data[index_cube];
                }
            }
        }
    } else {
        printf("read A file error!\n");
    }
    free(data);
    close(input_fd);
}

void init_x(char *filename_x, Data_Info info){
    int input_fd = open(filename_x, O_RDONLY);
    ll size = info.nx * info.ny * info.nz;
    double *data = (double *)malloc(size * sizeof(double));
    ll file_offset = info.x_start * info.tY * info.nz + info.y_start * info.nz + info.z_start;
    ll result = lseek(input_fd, file_offset * sizeof(double), SEEK_SET);
    if ((read(input_fd, data, sizeof(double) * size)) > 0){
        int i, j, k;
        ll size_cube = info.ny * info.nz;
        for (i = 0; i < info.nx; i++) {
            for (j = 0; j < info.ny; j++) {
                for (k = 0; k < info.nz; k++) {
                    ll index_cube = i * size_cube + j * info.nz + k; //row number of cube
                    //fill into x_0 in x y z order
                    x[index_cube] = data[index_cube];
                }
            }
        }

    } else {
        printf("read A file error!\n");
    }
    free(data);
    close(input_fd);
}

void A_m_vector(double * vector, int mrow, int mcol, double * result){
    int i, j;
    #pragma omp parallel private(i, j)
    {
        #pragma omp for schedule(static)
        for (i = 0; i < mrow; i++){
            int j_start = A_ptr[i];
            int j_end = A_ptr[i + 1];
            double tmp = 0;
            for (j = j_start; j < j_end; j++){
                tmp += A_value[j] * vector[A_col[j]];
            }
            result[i] = tmp;
        }
    }

}

void ILU_0(double *M, int rows){
    //只做方阵
    double j_row_V[30];
    int j_row_col[30];
    int j_count;
    int i;
    for (i = 0; i < rows; i++){
        int tmp_pos;
        int start_pos = A_ptr[i];
        int end_pos = A_ptr[i + 1];
        for (tmp_pos = start_pos; tmp_pos < end_pos; tmp_pos++){
            int j = A_col[tmp_pos];
            if (j > i - 1){
                break;
            }
            double a_jj;
            //-------
            j_count = 0;
            int start_pos_j = A_ptr[j];
            int end_pos_j = A_ptr[j + 1];
            int tmp_a_j;
            for (tmp_a_j = start_pos_j; tmp_a_j < end_pos_j; tmp_a_j++){
                int tmp_j_col = A_col[tmp_a_j];
                double val_j = M[tmp_a_j];
                if(tmp_j_col >= rows){
                    break;
                }
                j_row_V[j_count] = val_j;
                j_row_col[j_count] = tmp_j_col;
                j_count++;
                if(tmp_j_col == j){
                    a_jj = val_j;
                }
            }
            double a_ij = M[tmp_pos];
            a_ij = a_ij / a_jj;
            M[tmp_pos] = a_ij;
            int tmp_jpos;
            int j_row_search = 0;
            for(tmp_jpos = tmp_pos + 1; tmp_jpos < end_pos; tmp_jpos++){
                int k = A_col[tmp_jpos];
                if (k < rows && j_row_search < j_count){
                    //只使用方阵
                    while(j_row_search < j_count){
                        if(k == j_row_col[j_row_search]){
                            double a_jk = j_row_V[j_row_search];
                            double a_ik = M[tmp_jpos];
                            a_ik = a_ik - a_ij * a_jk;
                            M[tmp_jpos] = a_ik;
                            j_row_search++;
                            break;
                        } else {
                            if (k > j_row_col[j_row_search]){
                                j_row_search++;
                            } else {
                                break;
                            }
                        }
                    }
                } else {
                    break;
                }
            }
        }
    }
}

void ilu_solver(double *R, double *M, double *R_hat, int rows){
    int i;
    // 第一轮
    for (i = 0; i < rows; i++){
        int start_j = A_ptr[i];
        int end_j = A_ptr[i + 1];
        int j;
        double tmp = R[i];
        for (j = start_j; j < end_j; j++){
            int col_j = A_col[j];
            if (col_j < i){
                tmp -= M[j] * R_hat[col_j];
            } else {
                break;
            }
        }
        R_hat[i] = tmp;
    }
    // 第二轮
    for (i = rows - 1; i > -1; i--){
        int start_j = A_ptr[i];
        int end_j = A_ptr[i + 1];
        int j;
        double u_ii;
        double tmp = R_hat[i];
        for (j = start_j; j < end_j; j++){
            int col_j = A_col[j];
            if (col_j > i && col_j < rows){
                tmp -= M[j] * R_hat[col_j];
            } else {
                if (col_j == i){
                    u_ii = M[j];
                }
            }
        }
        R_hat[i] = tmp / u_ii;
    }
}

double avx_dot(double *A, double *B, int N) {
    vec temp = {0};
    N >>=2;
    vec *Av = (vec *)A;
    vec *Bv = (vec *)B;
    int i;
    for(i = 0; i < N; ++i) {
        temp += *Av * *Bv;
        Av++;
        Bv++;
    }
    union {
      vec tempv;
      double tempf[4];
    } u;

    u.tempv = temp;

    double dot = 0;
    for(i = 0; i < VECTOR_SIZE; ++i) {
        dot += u.tempf[i];
    }
    return dot;
}

void update_x(double *p_i, double alpha, int len){
    int i;
    #pragma omp parallel private(i)
    {
        #pragma omp for schedule(static)
        for (i = 0; i < len; i++){
            x[i] = x[i] + alpha * p_i[i];
        }
    }

}

int check_x(int rows){
    //需要进行 局部通讯以及 allreduce
    int i, j;
    double sum = 0;

    for (i = 0; i < rows; i++){
        int start_j = A_ptr[i];
        int end_j = A_ptr[i + 1];
        double b_i = b[i];
        double Ax_i = 0;
        for (j = start_j; j < end_j; j++){
            Ax_i += A_value[j] * x[A_col[j]];
        }
        sum += (b_i - Ax_i) * (b_i - Ax_i);
    }

    double g_sum;
    MPI_Allreduce(&sum, &g_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    g_sum = sqrt(g_sum);
    if (g_sum < RST){
        return 1;
    } else {
        return 0;
    }
}

void update_R(double *R, double alpha, double *Ap_i, int len){
    int i;
    #pragma omp parallel private(i)
    {
        #pragma omp for schedule(static)
        for (i = 0; i < len; i++){
            R[i] = R[i] - alpha * Ap_i[i];
        }
    }

}

void send_data(double *vector, int nx, int ny, int nz, int myrank){
    int i, j, k;
    //send side 0
    int size_yz = ny * nz;
    int start_x[8] = {0, 0, 0, 0, nx - 1, nx - 1, nx - 1, 0};
    int end_x[8] =   {1, 1, 1, nx, nx, nx, nx, nx};
    int start_y[8] = {0, 0, ny - 1, ny - 1, ny - 1, 0, 0, 0};
    int end_y[8] =   {1, ny, ny, ny, ny, ny, 1, 1};
    int side;

    #pragma omp parallel private(i, j, k, side)
    {
        #pragma omp for schedule(static)
        for (side = 0; side < 8; side++){
            if (myrank != node_rank[side]){
                // 需要发送
                int count = 0;
                double *start_index = send_buffer + buffer_offset[side];
                for (i = start_x[side]; i < end_x[side]; i++){
                    for (j = start_y[side]; j < end_y[side]; j++){
                        for (k = 0; k < nz; k++){
                            int index = i * size_yz + j * nz + k;
                            start_index[count] = vector[index];
                            count++;
                        }
                    }
                }
            }
        }
    }


    for (side = 0; side < 8; side++){
        // send data
        double *start_index = send_buffer + buffer_offset[side];
        MPI_Isend(start_index, node_dsize[side], MPI_DOUBLE, node_rank[side], target_side[side], MPI_COMM_WORLD, &send_request[side]);
    }

}

void recv_data(double *vector, int myrank){
    int side;
    for (side = 0; side < 8; side++){
        if (node_rank[side] != myrank){
            double *receive_buffer = vector + buffer_offset[side];
            MPI_Irecv(receive_buffer, node_dsize[side], MPI_DOUBLE, node_rank[side], side, MPI_COMM_WORLD, &recv_request[side]);
            //printf("recv ------- side : %d, node_dsize : %d, side : %d \n", side, node_dsize[side], side);
        }
    }
}

void wait_recv(int myrank){
    int side;
    MPI_Status status;
    for (side = 0; side < 8; side++){
        if (node_rank[side] != myrank){
            MPI_Wait(&recv_request[side], &status);
        }
    }
}

void gcr(Data_Info info){
    //init data
    double *R = (double *) malloc(info.len_x * sizeof(double));
    double *R_hat = (double *) malloc(info.len_x * sizeof(double)); // receive buffer is behind local values
    double *Ap = (double *) malloc(k * info.len_vb * sizeof(double));
    double *p = (double *) malloc(k * info.len_vb * sizeof(double));
    double *M = (double *) malloc(info.len_vb * 19 * sizeof(double));
    memcpy(M, A_value, info.len_vb * 19 * sizeof(double));
    //中间向量
    double *mid_result = (double *) malloc(info.len_vb * sizeof(double));
    double * recv_start;
    //init R
    //ssend and receive x
    send_data(x, info.nx, info.ny, info.nz, info.myrank);
    recv_start = x + info.len_vb;
    recv_data(recv_start, info.myrank);
    //waite recv done
    wait_recv(info.myrank);

    A_m_vector(x, info.len_vb, info.len_x, R); //ok

    int i;
    for (i = 0; i < info.len_vb; i++){
        R[i] = b[i] - R[i];
    }
    //init M
    ILU_0(M, info.len_vb);
    // get r hat
    ilu_solver(R, M, R_hat, info.len_vb);

    // init p_0
    memcpy(p, R_hat, info.len_vb * sizeof(double));
    // init Ap,  A * p_0 = A * R_hat , need communicate -------------------
    send_data(R_hat, info.nx, info.ny, info.nz, info.myrank);
    recv_start = R_hat + info.len_vb;
    recv_data(recv_start, info.myrank);
    wait_recv(info.myrank);

    A_m_vector(R_hat, info.len_vb, info.len_x, Ap);
    //printf("init Ap done\n");
    int steps = 0; //迭代的次数
    //开始迭代
    double *Ap_i;
    double *p_i;
    double beta[10];// k最大为10
    double Ap_idot[10];// k最大为10
    //printf("begin to iterate\n");
    Ap_i = Ap;
    p_i = p;
    while(1){

        //计算alpha 需要allreduce, R 需要局部通讯-------------------------------
        double l_numerator = avx_dot(R, Ap_i, info.len_vb);
        double l_denominator = avx_dot(Ap_i, Ap_i, info.len_vb);
        double denominator;
        MPI_Allreduce(&l_denominator, &denominator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double numerator;
        MPI_Allreduce(&l_numerator, &numerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        Ap_idot[steps % k] = denominator;
        double alpha = numerator / denominator;
        //开始更新 x, 不需要进行通讯
        update_x(p_i, alpha, info.len_vb);
        //开始check需要进行通讯---------------------
        send_data(x, info.nx, info.ny, info.nz, info.myrank);
        recv_start = x + info.len_vb;
        recv_data(recv_start, info.myrank);
        //waite recv done
        wait_recv(info.myrank);
        int judge = check_x(info.len_vb);
        //printf("check x done\n");
        if (judge == 1){
            printf("step: %d, x get\n", steps);
            break;
        } else {
            if (steps == 30){
                printf("failed\n");
                break;
            }
        }
        //开始更新R
        update_R(R, alpha, Ap_i, info.len_vb);
        //跟新 R_hat 不需要通讯
        ilu_solver(R, M, R_hat, info.len_vb);
        //R_hat = R;
        // 开始计算beta
        int sub_j;
        //计算 A * R_hat, 需要局部通讯到R_hat-----------------
        send_data(R_hat, info.nx, info.ny, info.nz, info.myrank);
        recv_start = R_hat + info.len_vb;
        recv_data(recv_start, info.myrank);
        wait_recv(info.myrank);
        A_m_vector(R_hat, info.len_vb, info.len_x, mid_result);

        for (sub_j = (steps / k) * k; sub_j <= steps; sub_j++){
            //计算分子的点积, 需要allreduce-------------------------
            double *Ap_j = Ap + (sub_j % k) * info.len_vb;
            l_numerator = avx_dot(mid_result, Ap_j, info.len_vb);
            MPI_Allreduce(&l_numerator, &numerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            denominator = Ap_idot[sub_j % k];
            beta[sub_j % k] = -numerator / denominator;
        }
        // 开始更新p_i
        steps++;
        p_i = p + (steps % k) * info.len_vb;
        int sub_i;

        for (sub_i = 0; sub_i < info.len_vb; sub_i++){
            int last_step = steps - 1;
            double tmp_result = 0;
            tmp_result += R_hat[sub_i];
            for (sub_j = (last_step / k) * k; sub_j <= last_step; sub_j++){
                double *tmp_ptr;
                tmp_ptr = p + (sub_j % k) * info.len_vb;
                tmp_result += beta[sub_j % k] * tmp_ptr[sub_i];
            }
            p_i[sub_i] = tmp_result;
        }


        //printf("updtae p done\n");
        //开始更新 Ap_i
        Ap_i = Ap + (steps % k) * info.len_vb;
        int last_step = steps - 1;
        for (sub_i = 0; sub_i < info.len_vb; sub_i++){
            double tmp_result = 0;
            tmp_result += mid_result[sub_i];
            for (sub_j = (last_step / k) * k; sub_j <= last_step; sub_j++){
                double *tmp_Ap;
                tmp_Ap = Ap + (sub_j % k) * info.len_vb;
                tmp_result += beta[sub_j % k] * tmp_Ap[sub_i];
            }
            Ap_i[sub_i] = tmp_result;
        }


        //printf("update Ap_i done\n");

    }
    free(R);
    free(R_hat);
    free(Ap);
    free(p);
    free(M);
    free(mid_result);
}

int main(int argc, char **argv) {
    int nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int PX = atoi(argv[1]);
    int PY = atoi(argv[2]);
    int PZ = atoi(argv[3]);
    int NX = atoi(argv[4]);;
    int NY = atoi(argv[5]);;
    int NZ = atoi(argv[6]);;
    char *filename_A = argv[7];
    char *filename_x0 = argv[8];
    char *filename_b = argv[9];
    if (PZ != 1){
        return 0;
    }
    struct timeval t1, t2;
    struct timeval t3, t4;
    Data_Info info = init_info(NX, NY, NZ, PX, PY, PZ);
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t3, NULL);
    init_A(filename_A, info);
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t4, NULL);
    if (info.myrank == 0) {
        printf("file A time  %.6lf\n", TIME(t3, t4));
    }
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t3, NULL);
    init_b(filename_b, info);
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t4, NULL);
    if (info.myrank == 0) {
        printf("file b time  %.6lf\n", TIME(t3, t4));
    }
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t3, NULL);
    init_x(filename_x0, info);
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t4, NULL);
    if (info.myrank == 0) {
        printf("file x_0 time  %.6lf\n", TIME(t3, t4));
    }
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t1, NULL);
    gcr(info);
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t2, NULL);
    if (info.myrank == 0) {
        printf("gcr time: %.6lf\n", TIME(t1, t2));
    }
    free(A_value);
    free(A_col);
    free(A_ptr);
    free(b);
    free(x);
    free(send_buffer);
    MPI_Finalize();
    return 0;
}