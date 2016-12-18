#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/fcntl.h>
#include <string.h>

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
//MPI_Request recv_request[8]; // receive request to check needed data has been received. not all are used
//MPI_Request send_request[8]; // send requets
int buffer_offset[8]; // the start index of par i to send
int total_nodes;
int node_dsize[8]; // this node communicates with 8 nodes at most, this array save the date size to transfer
int node_rank[8]; // save the node rank to communicate
int node_offset[8];// record the offset data to receive and send

int x_cor[19] = {0, -1, +1, 0, 0, +1, +1, -1, -1, 0, -1, +1, 0, 0, 0, -1, +1, 0, 0};
int y_cor[19] = {0, 0, 0, -1, +1, +1, -1, -1, +1, 0, 0, 0, -1, +1, 0, 0, 0, -1, +1};
int z_cor[19] = {0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, +1, +1, +1, +1, +1};

Data_Info init_info(int NX, int NY, int NZ, int PX, int PY, int PZ){
    Data_Info info;
    if (PX == 1, PY == 1, PZ == 1){
        //info.myrank
        info.nx = NX;
        info.ny = NY;
        info.nz = NZ;
        info.x_par = PX;
        info.y_par = PY;
        info.z_par = PZ;
        info.x_start = 0;
        info.x_end = NX;
        info.y_start = 0;
        info.y_end = NY;
        info.z_start = 0;
        info.z_end = NZ;
        info.tX = NX;
        info.tY = NY;
        info.tZ = NZ;
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
        x = (double *) malloc(info.len_x * sizeof(double));
        memset(x, 0, info.len_x * sizeof(double));
        send_buffer = (double *) malloc(size_buffer * sizeof(double));
        //init node info
        int i;
        for (i = 0; i < 8; i++){
            node_dsize[i] = -1;
            node_rank[i] = -1;
        }
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
            g_i += (info.tY >> 1);
        }
    }
    // x 循环回来
    g_i = (g_i + info.nx) % info.nx;
    if (g_i >= info.x_start && g_i < info.x_end && g_j >= info.y_start && g_j < info.y_end && g_k >= info.z_start && g_k < info.z_end){
        // the node is still in cube
        int l_i = g_i - info.x_start;
        int l_j = g_j - info.y_start;
        int l_k = g_k - info.z_start;
        return (l_i * (info.ny * info.nz) + l_j * info.nz + l_k);
    } else {
        // the node is outof cube
        int size = info.nx * info.ny * info.nz;
        if (i == -1){
            return (size + (j + 1) * info.nz + k);
        } else {
            if (i == info.nx){
                return (size + (info.ny + 2) * info.nz + info.nx * info.nz * 2 + (j + 1) * info.nz + k);
            } else {
                if (j == -1){
                    return (size + (info.ny + 2) * info.nz + info.nz * i + k);
                } else {
                    return (size + (info.ny + 2) * info.nz + info.nz * info.nx+ info.nz * i + k);
                }
            }
        }
    }
}

int pair_comparitor (const void* lhs, const void* rhs)
{
    return (((pair *)lhs)->col - ((pair *)rhs)->col);
}

void init_A(char *filename_A, Data_Info info){
    printf("begin to init A!\n");
    A_ptr[0] = 0;
    int input_fd = open(filename_A, O_RDONLY);
    ll size = info.nx * info.ny * info.nz;
    around_point *data = (around_point *)malloc(size * sizeof(around_point));
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
                        //printf("value : %.6lf , col : %d \n", A_value[count], A_col[count]);
                        count++;
                    }

                    A_ptr[index_cube + 1] = A_ptr[index_cube] + pair_count;
                }
            }
        }
    } else {
        printf("read A file error!\n");
    }

    close(input_fd);
}

void init_b(char *filename_b, Data_Info info){
    int input_fd = open(filename_b, O_RDONLY);
    ll size = info.nx * info.ny * info.nz;
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
    ll size = info.nx * info.ny * info.nz;
    v_struct *data = (v_struct *)malloc(size * sizeof(v_struct));
    if ((read(input_fd,data,sizeof(around_point)*size)) > 0){
        int i, j, k;
        ll size_cube = info.ny * info.nz;
        ll size_all = info.tY * info.tZ;
        for (i = 0; i < info.nx; i++){
            for (j = 0; j < info.ny; j++){
                for (k = 0; k < info.nz; k++){
                    ll index_cube = i * size_cube + j * info.nz + k; //row number of cube
                    ll index_all = (i + info.x_start) * size_all + (j + info.y_start) * info.tZ + k + info.z_start;
                    //fill into x_0 in x y z order
                    x[index_cube] = data[index_all].val;
                }
            }
        }
        // begin to deal side
        int count = info.nx * info.ny * info.nz;
        //line 0
        if (node_dsize[0] != -1){

        }
        if (node_dsize[1] != -1){

        }
        if (node_dsize[2] != -1){

        }
        if (node_dsize[3] != -1){

        }
        if (node_dsize[4] != -1){

        }
        if (node_dsize[5] != -1){

        }
        if (node_dsize[6] != -1){

        }
        if (node_dsize[7] != -1){

        }
        if (node_dsize[8] != -1){

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
            if (A_col[j] >= mcol){
                break;
            }
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
            if (i == 116484 && j == 130125){
                printf("get value %.10lf \n", matrix[index]);
            }
            return matrix[index];
        } else {
            if (A_col[index] > j){
                return 0;
            }
        }
    }
    return 0;
}

void ILU_0(double *M, int rows){
    //只做方阵
    int i;
    for (i = 0; i < rows; i++){
        int tmp_pos;
        int start_pos = A_ptr[i];
        int end_pos = A_ptr[i + 1];
        for (tmp_pos = start_pos; tmp_pos < end_pos; tmp_pos++){
            int k = A_col[tmp_pos];
            if (k >= i - 1){
                break;
            }
            double a_kk = get_value(k, k, M);
            //printf("a_{%d, %d} : %.12lf \n", k, k, a_kk);
            double a_ik = M[tmp_pos];
            /*
            if (i == 123323 && k == 130125){
                printf("pos : %d \n", tmp_pos);
                printf("a_{%d, %d} : %.12lf  A: %.12lf\n", i, k, M[tmp_pos], A_value[tmp_pos]);
            }
            */
            a_ik = a_ik / a_kk;
            M[tmp_pos] = a_ik;
            int tmp_jpos;
            for(tmp_jpos = tmp_pos + 1; tmp_jpos < end_pos; tmp_jpos++){
                int j = A_col[tmp_jpos];
                if (j >= rows){
                    //只使用方阵
                    break;
                }
                double a_ij = M[tmp_jpos];

                a_ij = a_ij - a_ik * get_value(k, j, M);
                M[tmp_jpos] = a_ij;

            }
        }
    }
}

void ilu_solver(double *R, double *M, double *R_hat, int rows){
    int i;
    // 第一轮
    R_hat[0] = R[0];
    for (i = 1; i < rows; i++){
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
        double u_ii = get_value(i, i, M);
        double tmp = R_hat[i];
        for (j = start_j; j < end_j; j++){
            int col_j = A_col[j];
            if (col_j > i && col_j < rows){
                tmp -= M[j] * R_hat[col_j];
            }
        }
        R_hat[i] = tmp / u_ii;
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
    sum = sqrt(sum);
    printf("sum : %.10lf \n", sum);
    if (sum < 0.000018789184606079){
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
    tmp = sqrt(tmp);

    printf("R sum : %.10lf \n", tmp);
}

void check_P(int rows, double *P){
    int i;
    double tmp = 0;
    for (i = 0; i < rows ; i++){
        tmp += P[i] * P[i];
    }
    tmp = sqrt(tmp);
    printf("P sum : %.10lf \n", tmp);
}

void update_R(double *R, double alpha, double *Ap_i, int len){
    int i;
    for (i = 0; i < len; i++){
        R[i] = R[i] - alpha * Ap_i[i];
    }
}

void gcr(Data_Info info, int k){
    printf("begin to cal gcr\n");
    //init data
    double *R = (double *) malloc(info.len_vb * sizeof(double));
    double *R_hat = (double *) malloc(info.len_x * sizeof(double)); // receive buffer is behind local values
    memset(R_hat, 0, info.len_x * sizeof(double));
    double *Ap = (double *) malloc(k * info.len_vb * sizeof(double));
    double *p = (double *) malloc(k * info.len_vb * sizeof(double));
    double *M = (double *) malloc(info.len_vb * 19 * sizeof(double));
    memcpy(M, A_value, info.len_vb * 19 * sizeof(double));
    //中间向量
    double *mid_result = (double *) malloc(info.len_vb * sizeof(double));
    //init R
    A_m_vector(x, info.len_vb, info.len_x, R); //ok
    int i;
    for (i = 0; i < info.len_vb; i++){
        R[i] = b[i] - R[i];
    }
    printf("init R done\n");
    //init M

    int tmp_p;
    for (tmp_p = 0; tmp_p < info.len_vb; tmp_p++){
        int s = A_ptr[tmp_p];
        int e = A_ptr[tmp_p + 1];
        int j;
        //printf("row : %d \n", tmp_p);
        for (j = s; j < e; j++){

            if (A_col[j] == 130125 && tmp_p == 130162){
                printf("A_%d%d : %.20lf\n", tmp_p, A_col[j], A_value[j]);
                printf("M_%d%d : %.20lf\n", tmp_p, A_col[j], M[j]);
                //printf("col: %d \n", tmp_p);
            }

            //printf("row: %d, col : %d \n", tmp_p, A_col[tmp_p]);
        }
    }
    ILU_0(M, info.len_vb);
    printf("init M done\n");
    // get r hat
    ilu_solver(R, M, R_hat, info.len_vb);
    printf("init R hat done\n");
    // init p_0
    memcpy(p, R_hat, info.len_vb * sizeof(double));
    // init Ap,  A * p_0 = A * R_hat , need communicate -------------------
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
        double numerator = vector_dot(R, Ap_i, info.len_vb);
        printf(" alpha numerator: %.30lf\n", numerator);
        double denominator = vector_dot(Ap_i, Ap_i, info.len_vb);
        printf(" alpha denominator: %.10lf\n", denominator);
        Ap_idot[steps % k] = denominator;
        double alpha = numerator / denominator;
        printf("cal alpha done, alpha: %.10lf\n", alpha);
        //开始更新 x, 不需要进行通讯
        update_x(p_i, alpha, info.len_vb);
        //printf("update x done\n");
        //开始check需要进行通讯---------------------
        int judge = check_x(info.len_vb);
        //printf("check x done\n");
        if (judge == 1){
            break;
        } else {
            printf("step: %d, not converge\n", steps);
        }
        //开始更新R
        update_R(R, alpha, Ap_i, info.len_vb);
        //check_R(info.len_vb, R);
        //printf("update R done\n");
        //跟新 R_hat 不需要通讯
        ilu_solver(R, M, R_hat, info.len_vb);
        //R_hat = R;
        //printf("update R hat done\n");
        // 开始计算beta
        int sub_j;
        //计算 A * R_hat, 需要局部通讯到R_hat-----------------
        A_m_vector(R_hat, info.len_vb, info.len_x, mid_result);
        for (sub_j = (steps / k) * k; sub_j <= steps; sub_j++){
            //计算分子的点积, 需要allreduce-------------------------
            double *Ap_j = Ap + (sub_j % k) * info.len_vb;
            numerator = vector_dot(mid_result, Ap_j, info.len_vb);
            denominator = Ap_idot[sub_j % k];
            beta[sub_j % k] = -numerator / denominator;
        }
        /*
        for (sub_j = 0; sub_j < k; sub_j++){
            printf("beta[%d] : %0.19lf; ", sub_j, beta[sub_j]);
        }
        */
        printf("\n");
        //printf("get beta done\n");
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
    }

}

int main(int argc, char **argv) {
    //MPI_Init(&argc, &argv);
    /*
    int NX = atoi(argv[1]);
    int NY = atoi(argv[2]);
    int NZ = atoi(argv[3]);
    int PX = atoi(argv[4]);
    int PY = atoi(argv[5]);
    int PZ = atoi(argv[6]);
    char *filename_A = argv[7];
    char *filename_x0 = argv[8];
    char *filename_b = argv[9];
     */
    int NX = 360;
    int NY = 180;
    int NZ = 38;
    int PX = 1;
    int PY = 1;
    int PZ = 1;
    char *filename_A = "./data/case_1bin/data_A.bin";
    char *filename_x0 = "./data/case_1bin/data_x0.bin";
    char *filename_b = "./data/case_1bin/data_b.bin";
    if (PZ != 1){
        return 0;
    }
    total_nodes = PX * PY * PZ;
    Data_Info info = init_info(NX, NY, NZ, PX, PY, PZ);
    printf("init info done!\n");
    init_A(filename_A, info);
    printf("init A done!\n");
    init_b(filename_b, info);
    printf("init b done!\n");
    init_x(filename_x0, info);
    printf("init x done!\n");
    gcr(info, 5);
    return 0;
}