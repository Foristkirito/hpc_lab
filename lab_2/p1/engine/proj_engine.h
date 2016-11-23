//
// Created by 夏鑫 on 22/11/2016.
//

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define MNRATIO 250
#define MAXP
#define NODES 16
typedef long long ll;
double alpha = 0.001; // 用来控制访存和网络的比率

//返回每一条边划分得到的份数, 1 表示 不划分
void get_partition(int nx, int ny, int nz, int *x_par, int *y_par, int *z_par){
    ll size_xy = nx * ny;
    ll size_xz = nx * nz;
    ll size_yz = ny * nz;
    int tmp_xp;
    int tmp_yp;
    int tmp_zp;
    int node_numbers = NODES;
    double mem_baseline = (ll)MNRATIO * (ll)nx * (ll)ny * (ll)nz * alpha;
    //printf("memory access baseline : %lf \n", mem_baseline);
    for (; node_numbers > 0; node_numbers--){
        ll mini_network = -1;
        for (tmp_xp = 1; tmp_xp <= node_numbers; tmp_xp++) {
            if (node_numbers % tmp_xp == 0) {
                int remain_x = node_numbers / tmp_xp;
                for (tmp_yp = 1; tmp_yp <= remain_x; tmp_yp++) {
                    if (remain_x % tmp_yp == 0) {
                        tmp_zp = remain_x / tmp_yp;
                        // begin to check baselines
                        ll network = node_numbers * 2 * ((tmp_xp - 1) * size_yz + (tmp_yp - 1) * size_xz + (tmp_zp - 1) * size_xy);
                        if (network <= mem_baseline) {
                            if (mini_network > 0){
                                if (network < mini_network){
                                    *x_par = tmp_xp;
                                    *y_par = tmp_yp;
                                    *z_par = tmp_zp;
                                    mini_network = network;
                                }
                            } else {
                                *x_par = tmp_xp;
                                *y_par = tmp_yp;
                                *z_par = tmp_zp;
                                mini_network = network;
                            }
                        }
                    }
                }
            }
        }
        if (mini_network > 0){
            //printf("network communication : %d \n", mini_network);
            return;
        }
    }

}
