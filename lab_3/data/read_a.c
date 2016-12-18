//
// Created by 夏鑫 on 13/12/2016.
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/fcntl.h>

typedef struct{
  int i;
  int j;
  int k;
  double val[19];
} around_point;

typedef long long ll;

int main(){
    size_t NX = 360;
    size_t NY = 180;
    size_t NZ = 38;
    size_t size = NX * NY * NZ;
    char file_name[100] = "./case_1bin/data_A.bin";
    int input_fd = open(file_name, O_RDONLY);
    around_point *data = (around_point *)malloc(size * sizeof(around_point));
    if ((read (input_fd,data,sizeof(around_point)*size)) > 0){
        int count = 0;
        size_t index;
        for (index = 0; index < size; index++){
            count++;
            if (count == 10000){
                printf("i : %d, j : %d, k : %d, val : %.12lf \n", data[index].i, data[index].j, data[index].k, data[index].val[5]);
                count = 0;
            }
        }
    }
    close(input_fd);
    return 0;
}
