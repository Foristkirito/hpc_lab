#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/fcntl.h>
#include <string.h>
#include <mpi.h>

#define TIME(a, b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))
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
  int i;
  int j;
  int k;
  double val[19];
} around_point;

typedef struct{
  int col;
  double val;
} pair;

typedef struct{
  int i;
  int j;
  int k;
  double val;
} v_struct;

double *A_value; // partial A
int *A_col;
int *A_ptr;
double *b; // partial b
double *x; // partial x, recieve buffer is behind the local x
double *send_buffer;//store the data to send, 8 segments, not all are used
int *g_map;//map behind index to global index
MPI_Request recv_request[8]; // receive request to check needed data has been received. not all are used
MPI_Request send_request[8]; // send requets
int buffer_offset[8]; // the start index of par i to send
int total_nodes;
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

int get_col(int i, int j, int k, Data_Info info, int *g_index){
    *g_index = 0;
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
        *g_index = g_i * info.tY * info.tZ + g_j * info.tZ + g_k;
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
    ll size = info.tX * info.tY * info.tZ;
    around_point *data = (around_point *)malloc(size * sizeof(around_point));

    A_ptr[0] = 0;
    if ((read(input_fd,data,sizeof(around_point)*size)) > 0){
        int i, j, k;
        ll size_cube = info.ny * info.nz;
        ll size_all = info.tY * info.tZ;
        int count = 0;
        for (i = 0; i < info.nx; i++){
            for (j = 0; j < info.ny; j++){
                for (k = 0; k < info.nz; k++){
                    pair data_tmp[19];
                    int pair_count = 0;
                    ll index_cube = i * size_cube + j * info.nz + k; //row number of cube
                    ll index_all = (i + info.x_start) * size_all + (j + info.y_start) * info.tZ + k + info.z_start;
                    //get col
                    int sub_num;
                    for (sub_num = 0; sub_num < 19; sub_num++){
                        double val = data[index_all].val[sub_num];
                        int g_index_v;
                        int col = get_col(i + x_cor[sub_num], j + y_cor[sub_num], k + z_cor[sub_num], info, &g_index_v);

                        if (col >= 0){
                            data_tmp[pair_count].val = val;
                            data_tmp[pair_count].col = col;

                            if (col >= size_cube * info.nx){
                                g_map[col - size_cube * info.nx] = g_index_v;
                            }
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
                    if (info.myrank == 0){
                        if (index_cube == 1231198){
                            printf("\n");
                        }
                    }
                    A_ptr[index_cube + 1] = A_ptr[index_cube] + pair_count;
                }
            }
        }
    } else {
        printf("read A file error!\n");
    }
    //free(data);
    close(input_fd);
}

void printf_x(Data_Info info){
    int i, j, k;
    //send side 0
    int nx = info.nx;
    int ny = info.ny;
    int nz = info.nz;
    int myrank = info.myrank;
    int size_yz = ny * nz;
    int start_x[8] = {0, 0, 0, 0, nx - 1, nx - 1, nx - 1, 0};
    int end_x[8] =   {1, 1, 1, nx, nx, nx, nx, nx};
    int start_y[8] = {0, 0, ny - 1, ny - 1, ny - 1, 0, 0, 0};
    int end_y[8] =   {1, ny, ny, ny, ny, ny, 1, 1};
    int side;
    int size = info.nx * info.ny * info.nz;
    int count = 0;
    for (side = 0; side < 8; side++){
        // 需要发送
        printf("-------------side %d \n", side);
        for (i = start_x[side]; i < end_x[side]; i++){
            for (j = start_y[side]; j < end_y[side]; j++){
                for (k = 0; k < nz; k++){
                    double val = x[count + size];
                    printf("side %d, i : %d, j : %d, k : %d, col : %d, val : %.10lf \n", side, i, j, k, count + size, val);
                    count++;
                }
            }
        }
    }
}

void init_b(char *filename_b, Data_Info info){
    int input_fd = open(filename_b, O_RDONLY);
    ll size = info.tX * info.tY * info.tZ;
    v_struct *data = (v_struct *)malloc(size * sizeof(v_struct));
    if ((read(input_fd,data,sizeof(around_point)*size)) > 0){
        int i, j, k;
        ll size_cube = info.ny * info.nz;
        ll size_all = info.tY * info.tZ;
        int count = 0;
        for (i = 0; i < info.nx; i++){
            for (j = 0; j < info.ny; j++){
                for (k = 0; k < info.nz; k++){
                    ll index_cube = i * size_cube + j * info.nz + k; //row number of cube
                    ll index_all = (i + info.x_start) * size_all + (j + info.y_start) * info.tZ + k + info.z_start;
                    //fill into b in x y z order
                    b[index_cube] = data[index_all].val;
                }
            }
        }
    } else {
        printf("read A file error!\n");
    }
    close(input_fd);
}

void init_x(char *filename_x, Data_Info info){
    int input_fd = open(filename_x, O_RDONLY);
    ll size = info.tX * info.tY * info.tZ;
    v_struct *data = (v_struct *)malloc(size * sizeof(v_struct));
    if ((read(input_fd,data,sizeof(around_point)*size)) > 0){
        int i, j, k;
        ll size_cube = info.ny * info.nz;
        ll size_all = info.tY * info.tZ;
        for (i = 0; i < info.nx; i++) {
            for (j = 0; j < info.ny; j++) {
                for (k = 0; k < info.nz; k++) {
                    ll index_cube = i * size_cube + j * info.nz + k; //row number of cube
                    ll index_all = (i + info.x_start) * size_all + (j + info.y_start) * info.tZ + k + info.z_start;
                    //fill into x_0 in x y z order
                    x[index_cube] = data[index_all].val;
                }
            }
        }

    } else {
        printf("read A file error!\n");
    }
    close(input_fd);
}

void A_m_vector(double * vector, int mrow, int mcol, double * result){
    int i, j;
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

double get_value(int i, int j, double *matrix){
    int start_j = A_ptr[i];
    int end_j = A_ptr[i + 1];
    int index;
    for (index = start_j; index < end_j; index++){
        if (A_col[index] == j){
            return matrix[index];
        } else {
            if (A_col[index] > j){
                return 0;
            }
        }
    }
    return 0;
}

void DILU(double *M, int rows){
    //只做方阵
    int i;
    for (i = 0; i < rows; i++){
        M[i] = get_value(i, i, A_value);
    }
    for (i = 0; i < rows; i++){
        int tmp_pos;
        int start_pos = A_ptr[i];
        int end_pos = A_ptr[i + 1];
        M[i] = 1 / M[i];
        int pos;
        for (pos = start_pos; pos < end_pos; pos++){
            int j = A_col[pos];
            if (j > i && j < rows){
                double a_ij = get_value(i, j, A_value);
                if (fabs(a_ij - 0) > 0.00000000001){
                    double a_ji = get_value(j, i, A_value);
                    if (fabs(a_ji - 0) > 0.00000000001){
                        M[j] = M[j] - a_ji * M[i] * a_ij;
                    }
                }
            }
        }
    }
}

void dilu_solver(double *R, double *M, double *R_hat, int rows){
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
                tmp -= A_value[j] * R_hat[col_j];
            } else {
                break;
            }
        }
        R_hat[i] = tmp / M[i];
    }
    // 第二轮
    for (i = rows - 1; i > -1; i--){
        int start_j = A_ptr[i];
        int end_j = A_ptr[i + 1];
        int j;
        double tmp = R_hat[i];
        double sum_tmp = 0;
        for (j = start_j; j < end_j; j++){
            int col_j = A_col[j];
            if (col_j > i && col_j < rows){
                sum_tmp += A_value[j] * R_hat[col_j];
            }
        }
        R_hat[i] = tmp - sum_tmp / M[i];
    }
}


double vector_dot(double *vec1, double *vec2, int len){
    double result = 0;
    int i;
    for (i = 0; i < len; i++){
        result += vec1[i] * vec2[i];
    }
    return result;
}

void update_x(double *p_i, double alpha, int len){
    int i;
    for (i = 0; i < len; i++){
        x[i] = x[i] + alpha * p_i[i];
    }
}

int check_x(int s, int rows, int myrank, Data_Info info){
    //需要进行 局部通讯以及 allreduce
    int i, j;
    double sum = 0;

    for (i = s; i < rows; i++){
        int start_j = A_ptr[i];
        int end_j = A_ptr[i + 1];
        double b_i = b[i];
        //double b_i = 0;
        double Ax_i = 0;
        for (j = start_j; j < end_j; j++){
            Ax_i += A_value[j] * x[A_col[j]];
        }
        sum += (b_i - Ax_i) * (b_i - Ax_i);
    }

    double g_sum;
    if(myrank == 0){
        printf("sum : %.10lf \n", sum);
    }
    MPI_Allreduce(&sum, &g_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    g_sum = sqrt(g_sum);
    if(myrank == 0){
        printf("g_sum : %.10lf \n", g_sum);
    }
    if (g_sum < 0.000018789184606079){
        return 1;
    } else {
        return 0;
    }
}

void check_R(int rows, double *R){
    int i;
    double tmp = 0;
    for (i = 0; i < rows ; i++){
        tmp += R[i] * R[i];
    }
}

void check_P(int rows, double *P){
    int i;
    double tmp = 0;
    for (i = 0; i < rows ; i++){
        tmp += P[i] * P[i];
    }
    tmp = sqrt(tmp);
    //printf("P sum : %.10lf \n", tmp);
}

void update_R(double *R, double alpha, double *Ap_i, int len){
    int i;
    for (i = 0; i < len; i++){
        R[i] = R[i] - alpha * Ap_i[i];
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

    for (side = 0; side < 8; side++){
        if (myrank != node_rank[side]){
            // 需要发送
            double *start_index = send_buffer + buffer_offset[side];
            int count = 0;
            for (i = start_x[side]; i < end_x[side]; i++){
                for (j = start_y[side]; j < end_y[side]; j++){
                    for (k = 0; k < nz; k++){
                        int index = i * size_yz + j * nz + k;
                        start_index[count] = vector[index];
                        count++;
                    }
                }
            }
            // send data
            MPI_Isend(start_index, node_dsize[side], MPI_DOUBLE, node_rank[side], target_side[side], MPI_COMM_WORLD, &send_request[side]);
            //printf("send ------ side : %d, node_dsize : %d, target side : %d\n", side, node_dsize[side], target_side[side]);
        }
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

void wait_send(int myrank){
    int side;
    MPI_Status status;
    for (side = 0; side < 8; side++){
        if (node_rank[side] != myrank){
            MPI_Wait(&send_request[side], &status);
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

void print_r_hat(double *v, int len){
    int i;
    for(i = 0; i < len; i++){
        printf("rows : %d, R_hat : %.20lf \n", i, v[i]);
    }
}

void check_sum_p(double *p, int len){
    int i;
    double tmp = 0;
    for (i = 0; i < len; i++){
        tmp += p[i] * p[i];
    }
    double g_sum;
    MPI_Allreduce(&tmp, &g_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    printf("---------------------check sum : %.10lf \n", g_sum);
}

void gcr(Data_Info info, int k){
    //printf("begin to cal gcr\n");
    //init data
    double *R = (double *) malloc(info.len_x * sizeof(double));
    if (R == NULL){
        printf("R null \n");
    }
    double *R_hat = (double *) malloc(info.len_x * sizeof(double)); // receive buffer is behind local values
    if (R_hat == NULL){
        printf("R_hat null \n");
    }
    memset(R_hat, 0, info.len_x * sizeof(double));
    double *Ap = (double *) malloc(k * info.len_vb * sizeof(double));
    if (Ap == NULL){
        printf("Ap null \n");
    }
    double *p = (double *) malloc(k * info.len_vb * sizeof(double));
    if (p == NULL){
        printf("p null \n");
    }
    double *M = (double *) malloc(info.len_vb * sizeof(double));
    if (M == NULL){
        printf("M null \n");
    }
    //memcpy(M, A_value, info.len_vb * 19 * sizeof(double));
    //中间向量
    double *mid_result = (double *) malloc(info.len_vb * sizeof(double));
    if (mid_result == NULL){
        printf("mid_result null \n");
    }
    double * recv_start;
    //init R
    //ssend and receive x
    send_data(x, info.nx, info.ny, info.nz, info.myrank);
    recv_start = x + info.len_vb;
    recv_data(recv_start, info.myrank);
    //waite recv done
    wait_recv(info.myrank);
    check_x(0, info.len_vb, info.myrank, info);

    A_m_vector(x, info.len_vb, info.len_x, R); //ok

    int i;
    for (i = 0; i < info.len_vb; i++){
        R[i] = b[i] - R[i];
    }
    //check_R(info.len_vb, R);
    //init M

    DILU(M, info.len_vb);

     */
    // get r hat
    dilu_solver(R, M, R_hat, info.len_vb);

    // init p_0
    memcpy(p, R_hat, info.len_vb * sizeof(double));
    // init Ap,  A * p_0 = A * R_hat , need communicate -------------------
    send_data(R_hat, info.nx, info.ny, info.nz, info.myrank);
    recv_start = R_hat + info.len_vb;
    recv_data(recv_start, info.myrank);
    wait_recv(info.myrank);

    A_m_vector(R_hat, info.len_vb, info.len_x, Ap);
    printf("init Ap done\n");
    int steps = 0; //迭代的次数
    //开始迭代
    double *Ap_i;
    double *p_i;
    double beta[10];// k最大为10
    double Ap_idot[10];// k最大为10
    printf("begin to iterate\n");
    Ap_i = Ap + (steps % k) * info.len_vb;
    p_i = p + (steps % k) * info.len_vb;
    while(1){

        //计算alpha 需要allreduce, R 需要局部通讯-------------------------------
        double l_numerator = vector_dot(R, Ap_i, info.len_vb);
        printf(" alpha numerator_l: %.30lf\n", l_numerator);
        double l_denominator = vector_dot(Ap_i, Ap_i, info.len_vb);
        printf(" alpha denominator local: %.10lf\n", l_denominator);
        double denominator;
        MPI_Allreduce(&l_denominator, &denominator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double numerator;
        MPI_Allreduce(&l_numerator, &numerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        printf(" alpha denominator: %.10lf\n", denominator);
        printf(" alpha numerator: %.10lf\n", numerator);
        Ap_idot[steps % k] = denominator;
        double alpha = numerator / denominator;
        printf("cal alpha done, alpha: %.10lf\n", alpha);
        //开始更新 x, 不需要进行通讯
        update_x(p_i, alpha, info.len_vb);
        //printf("update x done\n");
        //开始check需要进行通讯---------------------
        send_data(x, info.nx, info.ny, info.nz, info.myrank);
        recv_start = x + info.len_vb;
        recv_data(recv_start, info.myrank);
        //waite recv done
        wait_recv(info.myrank);
        int judge = check_x(0, info.len_vb, info.myrank, info);

        //printf("check x done\n");
        if (judge == 1){
            printf("step: %d, x get\n", steps);
            break;
        } else {
            printf("step: %d, not converge\n", steps);
        }
        //开始更新R
        update_R(R, alpha, Ap_i, info.len_vb);
        //跟新 R_hat 不需要通讯
        dilu_solver(R, M, R_hat, info.len_vb);
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
            l_numerator = vector_dot(mid_result, Ap_j, info.len_vb);
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
        for (sub_i = 0; sub_i < info.len_vb; sub_i++){
            int last_step = steps - 1;
            double tmp_result = 0;
            tmp_result += mid_result[sub_i];
            for (sub_j = (last_step / k) * k; sub_j <= last_step; sub_j++){
                double *tmp_Ap;
                tmp_Ap = Ap + (sub_j % k) * info.len_vb;
                tmp_result += beta[sub_j % k] * tmp_Ap[sub_i];
            }
            Ap_i[sub_i] = tmp_result;
        }
        printf("update Ap_i done\n");
        if (steps == 150){
            printf("failed \n");
            break;
        }
    }
}

int main(int argc, char **argv) {
    //MPI_Init(&argc, &argv);
    /*
    int NX = atoi(argv[1]);
    int NY = atoi(argv[2]);
    int NZ = atoi(argv[3]);

    char *filename_A = argv[7];
    char *filename_x0 = argv[8];
    char *filename_b = argv[9];
     */
    int nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int NX = 360;
    int NY = 180;
    int NZ = 38;
    int PX = atoi(argv[1]);
    int PY = atoi(argv[2]);
    int PZ = atoi(argv[3]);
    total_nodes = PX * PY * PZ;
    char *filename_A = "./data/case_1bin/data_A.bin";
    char *filename_x0 = "./data/case_1bin/data_x0.bin";
    char *filename_b = "./data/case_1bin/data_b.bin";
    if (PZ != 1){
        return 0;
    }
    total_nodes = PX * PY * PZ;
    Data_Info info = init_info(NX, NY, NZ, PX, PY, PZ);
    //printf("init info done!\n");
    init_A(filename_A, info);
    //printf("init A done!\n");
    init_b(filename_b, info);
    //printf("init b done!\n");
    init_x(filename_x0, info);
    //printf("init x done!\n");
    struct timeval t1, t2;
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t1, NULL);
    gcr(info, 5);
    MPI_Barrier(MPI_COMM_WORLD), gettimeofday(&t2, NULL);
    if (info.myrank == 0) {
        printf("Total time: %.6lf\n", TIME(t1, t2));
    }
    MPI_Finalize();
    return 0;
}