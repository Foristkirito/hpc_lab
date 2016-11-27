#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../proj_engine.h"

int main(){
    int nx;
    int nz;
    int ny;
    scanf("%d", &nx);
    scanf("%d", &ny);
    scanf("%d", &nz);
    int x_par;
    int y_par;
    int z_par;
    get_partition(nx, ny, nz, &x_par, &y_par, &z_par);
    printf("x_par %d, y_par %d, z_par %d \n", x_par, y_par, z_par);
    return 0;
}