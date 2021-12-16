#include<stdio.h>
#include"TrackingFunctions.h"
#include"EM_Tuning_Tracking.h"

int main(int argc, char **argv)
{
    if(argc!=3){
        fprintf(stdout,"\n./run pic mv0\n");
        return 1;
    }
    /*read mv0*/
    FILE *mv0in;
    mv0in = fopen(argv[2],"r");
    struct State_Trace_Origin tracein;
    struct State_Trace tracemv0;
    int l;
    fscanf(mv0in,"%d",&l);
    Malloc_State_Trace_Origin(&tracein, l, 1);
    tracein.IDnum = l;
    int i,j;
    double x;
    for(i=0;i<l;i++){
	fscanf(mv0in,"%lf",&x);
	tracein.Y[i] = x;
	fscanf(mv0in,"%lf",&x);
	tracein.X[i] = x;
	fscanf(mv0in,"%lf",&x);
	tracein.I[i] = x;
	tracein.T[i] = 0;
	tracein.ID[i] = i;
    }
    fscanf(mv0in,"%lf",&x);
    fscanf(mv0in,"%lf",&x);
    fscanf(mv0in,"%lf",&x);
    tracein.no[0] = x;
    fclose(mv0in);
    Make_State_Trace(&tracein, &tracemv0);
    /*read parameters*/
    int bsize;
    int psfdecay;
    double sigma;
    double D;
    double alpha;
    int itemax;
    fprintf(stdout,"\ninput boundary size\n");
    fscanf(stdin,"%d",&bsize);
    fprintf(stdout,"\ninput maximal psf width, out of this psf are zeros\n");
    fscanf(stdin,"%d",&psfdecay);
    fprintf(stdout,"\ninput psf width, gaussian sigma\n");
    fscanf(stdin,"%lf",&sigma);
    fprintf(stdout,"\ninput diffusion sigma\n");
    fscanf(stdin,"%lf",&D);
    fprintf(stdout,"\ninput intensity linking rate\n");
    fscanf(stdin,"%lf",&alpha);
    fprintf(stdout,"\ninput maximal number of iterations\n");
    fscanf(stdin,"%d",&itemax);
    /*read the picture*/
    FILE *picin;
    picin = fopen(argv[1],"r");
    struct Img_Series picinput;
    int s1,s2;
    fscanf(picin,"%d",&s1);
    fscanf(picin,"%d",&s2);
    picinput.s1 = s1;
    picinput.s2 = s2;
    picinput.bsize = bsize;
    picinput.Tmax = 1;
    picinput.b = (double**)malloc(sizeof(double*)*1);
    picinput.b[0] = (double*)malloc(sizeof(double)*s1*s2);
    int idx;
    for(i=0;i<s1;i++){
	for(j=0;j<s2;j++){
	    idx = i+j*s1;
	    fscanf(picin,"%lf",picinput.b[0]+idx);
	}
    }
    Make_ExtraBoundary(&picinput);
    fclose(picin);
    /*do the iterations*/
    struct State_Trace statethis;
    Malloc_State_Trace(&tracemv0,&statethis);
    Deepcopy_State_Trace(&tracemv0,&statethis);
    struct EM_Trace_Run runpara;
    /*assign exptable*/
    double *exptable;
    double *expx;
    exptable=(double *)malloc(sizeof(double)*Expmin*Expinterval);
    expx=(double *)malloc(sizeof(double)*Expmin*Expinterval);
    int k;
    double dx;
    k=Expmin*Expinterval;
    dx=1.0/Expinterval;
    for(i=0;i<k;i++){
	expx[i]=-dx*i;
	exptable[i]=exp(-dx*i);
    }
    /*assign runpara*/
    runpara.s1 = s1;
    runpara.s2 = s2;
    runpara.T = picinput.Tmax;
    runpara.sig = sigma;
    runpara.psfdecay = psfdecay;
    runpara.bsize = bsize;
    runpara.ITmax = itemax;
    runpara.D = D;
    runpara.alpha = alpha;
    runpara.lambda = 0.0;
    runpara.exptable = exptable;
    runpara.expx = expx;
    /*do iterations*/
    for(i=0;i<itemax;i++){
	EM_Tuning_Step(&runpara, &picinput, &statethis, false);
	/*output the intermedie data*/
    }
    for(i=0;i<statethis.N;i++){
	fprintf(stdout,"\n%lf %lf\n",statethis.X[i][0],statethis.Y[i][0]);
    }

    /*free the space*/
    Free_State_Trace_Origin(&tracein);
    Free_State_Trace(&tracemv0);
    Free_State_Trace(&statethis);
    Free_Img_Series(&picinput);
    free(exptable);
    free(expx);
    return 1;
}


