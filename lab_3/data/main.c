//
// Created by 夏鑫 on 13/12/2016.
//

//write data struct to binary file

#include <stdio.h>
#include <stdlib.h>

typedef struct{
  int i;
  int j;
  int k;
  double val[19];
} around_point;

typedef struct{
  int i;
  int j;
  int k;
  double val;
} b_struct;

typedef long long ll;

int main(){
    int NX = 360;
    int NY = 180;
    int NZ = 38;
    int size = NX * NY * NZ;
    /*
    char file_name[100] = "./case_360x180x38/A.txt";
    FILE *fp=fopen(file_name,"r");
    if(fp==NULL){
        printf("can not open file!\n");
        exit(0);
    }
    around_point *data = (around_point *)malloc(size * sizeof(around_point));
    int i,j,k,num;
    double val;
    int index = 0;
    int sub_index;
    int count = 0;
    printf("begin to deal \n");
    for (;index < size; index++){
        for (sub_index = 0; sub_index < 19; sub_index++){
            fscanf(fp,"%d",&i); fscanf(fp,"%d",&j);
            fscanf(fp,"%d",&k); fscanf(fp,"%d",&num);
            fscanf(fp,"%lf",&val);
            data[index].val[sub_index] = val;
            if (sub_index == 0){
                data[index].i = i;
                data[index].j = j;
                data[index].k = k;
                if (count == 10000){
                    printf("i: %d, j : %d, k : %d \n", i, j, k);
                    count = 0;
                }
            }

        }
        count++;
    }
    fclose(fp);
    FILE *fw=fopen("data_A.bin","wb");
    if(fw==NULL){
        printf("can't open the write file\n");
        exit(0);
    }
    fwrite(data, sizeof(around_point), size, fw);
    fclose(fw);
    free(data);
    */
    char file_name[100] = "./case_360x180x38/x0.txt";
    FILE *fp=fopen(file_name,"r");
    if(fp==NULL){
        printf("can not open file!\n");
        exit(0);
    }
    b_struct *data = (b_struct *)malloc(size * sizeof(b_struct));
    int i,j,k;
    double val;
    int index = 0;
    int count = 0;
    printf("begin to deal \n");
    for (;index < size; index++){
        fscanf(fp,"%d",&i);
        fscanf(fp,"%d",&j);
        fscanf(fp,"%d",&k);
        fscanf(fp,"%lf",&val);
        data[index].i = i;
        data[index].j = j;
        data[index].k = k;
        data[index].val = val;
        count++;
        if (count == 10000){
            printf("i : %d, j : %d, k : %d, val : %.12lf \n", i, j, k, val);
            count = 0;
        }
    }
    fclose(fp);
    FILE *fw=fopen("./case_1bin/data_x0.bin","wb");
    if(fw==NULL){
        printf("can't open the write file\n");
        exit(0);
    }
    fwrite(data, sizeof(b_struct), size, fw);
    fclose(fw);
    free(data);
    return 0;
}