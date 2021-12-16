#include<stdio.h>
#include"TrackingFunctions.h"

int main(int argc, char **argv)
{
    double *a;
    double *ans;
    a = (double*)malloc(sizeof(double)*10);
    ans = (double*)malloc(sizeof(double)*10);
    int i;
    for(i=0;i<10;i++){
	a[i] = 1.0*i;
    }
    for(i=0;i<100000000;i++){
	Intensity_Update_Homotopy(10, a, 1000, ans);
    }
    free(a);
    free(ans);
    return 1;
}


