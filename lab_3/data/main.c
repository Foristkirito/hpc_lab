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
    int NX = 720;
    int NY = 360;
    int NZ = 38;
    ll size = NX * NY * NZ;


    /*
    char file_name[100] = "./case_720x360x38/A.txt";
    FILE *fp=fopen(file_name,"r");
    if(fp==NULL){
        printf("can not open file!\n");
        exit(0);
    }
    double *data = (double *)malloc(size * 19 * sizeof(double));
    int i,j,k,num;
    double val;
    ll index = 0;
    ll sub_index;
    ll count = 0;
    printf("begin to deal \n");
    for (;index < size; index++){
        for (sub_index = 0; sub_index < 19; sub_index++){
            fscanf(fp,"%d",&i); fscanf(fp,"%d",&j);
            fscanf(fp,"%d",&k); fscanf(fp,"%d",&num);
            fscanf(fp,"%lf",&val);
            data[19 * index + sub_index] = val;
            if (sub_index == 0){
                if (count == 10000){
                    printf("i: %d, j : %d, k : %d \n", i, j, k);
                    count = 0;
                }
            }

        }
        count++;
    }
    fclose(fp);
    FILE *fw=fopen("./case_2bin/data_A_v1.bin","wb");
    if(fw==NULL){
        printf("can't open the write file\n");
        exit(0);
    }
    fwrite(data, sizeof(double), size * 19, fw);
    fclose(fw);
    free(data);
    */


    char file_name[100] = "./case_720x360x38/x0.txt";
    FILE *fp=fopen(file_name,"r");
    if(fp==NULL){
        printf("can not open file!\n");
        exit(0);
    }
    double *data = (double *)malloc(size * sizeof(double));
    int i,j,k;
    double val;
    ll index = 0;
    ll count = 0;
    printf("begin to deal \n");
    for (;index < size; index++){
        fscanf(fp,"%d",&i);
        fscanf(fp,"%d",&j);
        fscanf(fp,"%d",&k);
        fscanf(fp,"%lf",&val);
        data[index] = val;
        count++;
        if (count == 10000){
            printf("i : %d, j : %d, k : %d, val : %.12lf \n", i, j, k, val);
            count = 0;
        }
    }
    fclose(fp);
    FILE *fw=fopen("./case_2bin/data_x0_v1.bin","wb");
    if(fw==NULL){
        printf("can't open the write file\n");
        exit(0);
    }
    fwrite(data, sizeof(double), size, fw);
    fclose(fw);
    free(data);

    return 0;
}