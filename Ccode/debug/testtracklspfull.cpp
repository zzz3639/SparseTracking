#include<stdio.h>
#include"TrackingFunctions.h"
#include"EM_Tuning_Tracking.h"

void Distant_State_Trace(State_Trace *st1, State_Trace *st2, double *distI, double *distX)
{
    *distI = 0.0;
    *distX = 0.0;
    int i,j;
    for(i=0;i<st1->N;i++){
	for(j=0;j<st1->n[i];j++){
	    *distI = *distI + fabs( st1->I[i][j] - st2->I[i][j] );
	    *distX = *distX + fabs( st1->X[i][j] - st2->X[i][j] );
	    *distX = *distX + fabs( st1->Y[i][j] - st2->Y[i][j] );
	}
    }
    for(i=0;i<st1->L;i++){
	*distI = *distI + fabs( st1->no[i] - st2->no[i] );
    }
    return;
}

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
    int T;
    fscanf(mv0in,"%d",&l);
    fscanf(mv0in,"%d",&T);
    Malloc_State_Trace_Origin(&tracein, l, T);
    int i,j;
    double x;
    for(i=0;i<l;i++){
	fscanf(mv0in,"%lf",&x);
	tracein.Y[i] = x;
	fscanf(mv0in,"%lf",&x);
	tracein.X[i] = x;
	fscanf(mv0in,"%lf",&x);
	tracein.I[i] = x;
        fscanf(mv0in,"%d",tracein.T+i);
	fscanf(mv0in,"%d",tracein.ID+i);
    }
    for(i=0;i<T;i++){
	fscanf(mv0in,"%lf",&x);
        tracein.no[i] = x;
    }
    fclose(mv0in);
    Reassign_State_Trace_Origin(&tracein);
    Make_State_Trace(&tracein, &tracemv0);
    /*read parameters*/
    int bsize;
    int psfdecay;
    double sigma;
    double D;
    double alpha;
    double *lambda;
    double *deltaI;
    int itemax;
    double molfiltervalue;
    lambda = (double*)malloc(sizeof(double)*tracein.IDnum);
    deltaI = (double*)malloc(sizeof(double)*tracein.IDnum);
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
    fprintf(stdout,"\ninput molecule filter intensity\n");
    fscanf(stdin,"%lf",&molfiltervalue);
    fprintf(stdout,"\ninput sparse penalty lambda, %d doubles\n",tracein.IDnum);
    for(i=0;i<tracein.IDnum;i++){
	fscanf(stdin,"%lf",lambda+i);
    }
    fprintf(stdout,"\ninput sparse penalty deltaI, %d doubles\n",tracein.IDnum);
    for(i=0;i<tracein.IDnum;i++){
	fscanf(stdin,"%lf",deltaI+i);   
    }
    fprintf(stdout,"\n%d molecules initiated, noise %lf, %d iterations\n",tracemv0.N, tracein.no[0], itemax);
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
    fscanf(picin,"%d",&(picinput.Tmax));
    picinput.b = (double**)malloc(sizeof(double*)*picinput.Tmax);
    for(i=0;i<picinput.Tmax;i++){
	picinput.b[i] = (double*)malloc(sizeof(double)*s1*s2);
    }
    int idx;
    int t;
    for(t=0;t<picinput.Tmax;t++){
	for(i=0;i<s1;i++){
	    for(j=0;j<s2;j++){
		idx = i+j*s1;
		fscanf(picin,"%lf",picinput.b[t]+idx);
	    }
	}
    }
    Make_ExtraBoundary(&picinput);
    fclose(picin);
    for(i=0;i<tracemv0.N;i++){
	tracemv0.sp[i] = 's';
	tracemv0.lambda[i] = lambda[i];
        tracemv0.deltaI[i] = deltaI[i];
    }
    /*prepare the iterations*/
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
    runpara.exptable = exptable;
    runpara.expx = expx;
    /*do iterations*/
    struct State_Trace statebuf;
    Malloc_State_Trace(&statethis,&statebuf);
    Deepcopy_State_Trace(&statethis,&statebuf);
    double distI,distX;
    double llh;
    for(i=0;i<itemax;i++){
        /*output the intermedie data*/
        //if(i%1==0){
	//    Distant_State_Trace(&statethis,&statebuf,&distI,&distX);
        //    llh = EM_Tuning_LogLikelihood(&runpara, &picinput, &statethis, false);
        //    fprintf(stdout,"\niteration: %d, convergence indicator: %lf %lf, loglikelihood: %f\n",i,distI,distX,llh);
        //    Deepcopy_State_Trace(&statethis,&statebuf);
        //}
	EM_Tuning_Step(&runpara, &picinput, &statethis, false);
    }
    bool *moldie;
    moldie = (bool*)malloc(sizeof(bool)*statethis.N);
    for(i=0;i<statethis.N;i++){
	x = 0.0;
	for(j=0;j<statethis.n[i];j++){
	    if(x<statethis.I[i][j])
                x=statethis.I[i][j];
	}
	if(x<molfiltervalue){
	    moldie[i] = true;
	}
	else{
	    moldie[i] = false;
	}
    }
    k = 0;
    for(i=0;i<statethis.N;i++){
	if(moldie[i]){
	}
	else{
	    for(j=0;j<statethis.n[i];j++){
		fprintf(stdout,"\n%lf %lf %lf %d %d\n",statethis.X[i][j],statethis.Y[i][j],statethis.I[i][j],statethis.T[i][j],k);
	    }
	    k++;
	}
    }
    for(i=0;i<statethis.L;i++){
        fprintf(stdout,"\n%lf\n",statethis.no[i]);
    }
    fprintf(stdout,"\n%d\n",k);
    /*free the space*/
    Free_State_Trace_Origin(&tracein);
    Free_State_Trace(&tracemv0);
    Free_State_Trace(&statethis);
    Free_Img_Series(&picinput);
    free(exptable);
    free(expx);
    free(moldie);
    free(lambda);
    free(deltaI);
    return 1;
}


