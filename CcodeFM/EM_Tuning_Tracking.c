#include"EM_Tuning_Tracking.h"


/*sig: psf width; bsize: extra boundary size; psfdecay: how far away over which psf decay to zero*/ 
/*D guassian sigma of molecule diffusion velocity; alpha: log laplace parameter of intensity connections; lambda: sparse penalty*/
/*state0: initial molecule lists; */
/*ITmax: max iteration number; if after llhwindow iterations, log likelihood doesn't change more than llhzero, stop the iteration*/ 
/*fix_frame0: the positions and intensities of frame0 is fixed or not*/
/*exptable and expx are allocated and computed outside*/
/*return value: log likelihood*/
double EM_Tuning(struct Img_Series *mv, double sig, int psfdecay, int bsize, double D, double alpha, double lambda, struct State_Trace *state0, int ITmax, double llhzero, int llhwindow, bool fix_frame0, double *exptable, double *expx)
{
    int i,j,k;
    /*Initialize struct EM_Trace_Run*/
    struct EM_Trace_Run run;
    run.s1 = mv->s1;
    run.s2 = mv->s2;
    run.T = mv->Tmax;
    run.sig = sig;
    run.psfdecay = psfdecay;
    run.bsize = bsize;
    run.ITmax = ITmax;
    run.D = D;
    run.alpha = alpha;
    run.exptable = exptable;
    run.expx = expx;
    /*vector record log likelihood along iterations*/
    double *llh;
    llh = (double *)malloc(sizeof(double)*(ITmax+1));
    llh[0] = EM_Tuning_LogLikelihood(&run, mv, state0, fix_frame0);

    /*do EM iterations*/
    int t;
    /*deep copy state0 to statenow*/
    struct State_Trace statenow;
    Malloc_State_Trace(state0, &statenow);
    Deepcopy_State_Trace(state0, &statenow);
    for(i=0;i<ITmax;i++){
	k = i+1;
	llh[k] = EM_Tuning_Step(&run, mv, &statenow, fix_frame0);
	if(k>llhwindow&&fabs(llh[k]-llh[k-llhwindow])<llhzero){
	    break;
	}
    }
    /*free allocated spaces inside this function*/
    Free_State_Trace(&statenow);
    double ans;
    ans = llh[k];
    free(llh);
    return ans;
}


/*the sparse penalty is either L1 or LSP, L1: lambda*I, LSP: lambda*log(1+I/deltaI)*/
double EM_Tuning_Step(struct EM_Trace_Run *run, struct Img_Series *mv, struct State_Trace *statethis, bool fix_frame0)
{
/*Initialize*/
    int i,j,k;
    int s1;
    int s2;
    s1 = run->s1;
    s2 = run->s2;
    int sb1;
    int sb2;
    sb1 = s1 + run->bsize + run->bsize;
    sb2 = s2 + run->bsize + run->bsize;
    int sbs;
    sbs = sb1*sb2;
    /*save the images, the boundary value of bb well be changed, we don't want to mess up the value inside run->bb */
    double **bb;
    bb = (double **)malloc(sizeof(double*)*(mv->Tmax));
    for(i=0;i<mv->Tmax;i++){
	bb[i] = (double*)malloc(sizeof(double)*sbs);
    }
    for(i=0;i<mv->Tmax;i++){
	for(j=0;j<sbs;j++){
	    bb[i][j] = mv->bb[i][j];
	}
    }
    /*save psf of the molecules*/
    struct Sparse **W;
    W = (struct Sparse **)malloc(sizeof(struct Sparse *)*(statethis->N));
    for(i=0;i<statethis->N;i++){
	W[i] = (struct Sparse *)malloc(sizeof(struct Sparse)*(statethis->n[i]));
    }
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i];j++){
	    Sparse_malloc(W[i]+j,(2+run->psfdecay+run->psfdecay)*(2+run->psfdecay+run->psfdecay));
	    W[i][j].n = 0;
	}
    }
    /*save psfno*/
    double **Wn;
    Wn = (double **)malloc(sizeof(double*)*(mv->Tmax));
    for(i=0;i<mv->Tmax;i++){
	Wn[i] = (double *)malloc(sizeof(double)*sbs);
    }
    /*save the sum up of psfs, expected images*/
    double **sumpsf;
    sumpsf = (double **)malloc(sizeof(double*)*(mv->Tmax));
    for(i=0;i<mv->Tmax;i++){
	sumpsf[i] = (double *)malloc(sizeof(double)*sbs);
    }
    /*save itemized inverse of sumpsf*/
    double **sumpsfinv;
    sumpsfinv = (double **)malloc(sizeof(double*)*(mv->Tmax));
    for(i=0;i<mv->Tmax;i++){
	sumpsfinv[i] = (double *)malloc(sizeof(double)*sbs);
    }
    /*save sum of W and W.*X for each molecule */
    struct State_Trace statetemp;
    Malloc_State_Trace(statethis, &statetemp);
    /*to save the matrix for position update*/
    struct EM_Matrix *Inv_Matrix;
    Inv_Matrix = (struct EM_Matrix*)malloc(sizeof(struct EM_Matrix)*statethis->N);
    for(i=0;i<statethis->N;i++){
        Malloc_EM_Matrix(Inv_Matrix+i, statethis->n[i], statethis->n[i]);
    }

/*E-step, compute and sum up psfs*/
    int idx;
    int T;
    /*compute psfs*/
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i];j++){
	    PSFGauss_Native(W[i]+j, run->exptable, run->expx, 0, sb1-1, 0, sb2-1, statethis->X[i][j], statethis->Y[i][j], run->sig, run->psfdecay);
	    for(k=0;k<W[i][j].n;k++){
	        W[i][j].V[k] = W[i][j].V[k] * statethis->I[i][j];
	    }
	}
    }
    /*compute psfno*/
    for(i=0;i<mv->Tmax;i++){
	PSFno(Wn[i],sbs,statethis->no[i]);
    }
    /*sum up psfs*/
    for(i=0;i<mv->Tmax;i++){
	for(j=0;j<sbs;j++){
	    sumpsf[i][j] = Wn[i][j];
	}
    }
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i];j++){
	    T = statethis->T[i][j];
	    for(k=0;k<W[i][j].n;k++){
		idx = W[i][j].X[k] + sb1*W[i][j].Y[k];
		sumpsf[T][idx] += W[i][j].V[k];
	    }
	}
    }
    /*reconstruct boundary part of the images*/
    for(i=0;i<mv->Tmax;i++){
	for(j=0;j<sbs;j++){
	    if(mv->bound[i][j]){
		continue;
	    }
	    sumpsfinv[i][j] = 1.0/sumpsf[i][j];
	}
    }
    /*count expected number of photons each molecule contribute at each pixel*/
    for(i=0;i<mv->Tmax;i++){
	for(j=0;j<sbs;j++){
            if(mv->bound[i][j]){
                continue;
            }
	    Wn[i][j] = Wn[i][j] * sumpsfinv[i][j] * bb[i][j];
	}
    }
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i];j++){
	    T = statethis->T[i][j];
	    for(k=0;k<W[i][j].n;k++){
		idx = W[i][j].X[k] + sb1*W[i][j].Y[k];
                if(mv->bound[T][idx]){
                    W[i][j].V[k] = W[i][j].V[k] * (1.0 + statethis->lambda[i]);
                    continue;
                }
		W[i][j].V[k] = W[i][j].V[k] * sumpsfinv[T][idx] * bb[T][idx];
	    }
	}
    }
/*M-step, update tracking states*/
    /*update noise values*/
    for(i=0;i<mv->Tmax;i++){
        statethis->no[i] = 0.0;
        for(j=0;j<sbs;j++){
            statethis->no[i] += Wn[i][j];
        }
    }
    /* sum up W and W.*X */
    for(i=0;i<statethis->N;i++){
        for(j=0;j<statethis->n[i];j++){
            statetemp.I[i][j] = 0.0;
            statetemp.X[i][j] = 0.0;
            statetemp.Y[i][j] = 0.0;
            T = statethis->T[i][j];
            for(k=0;k<W[i][j].n;k++){
                statetemp.I[i][j] += W[i][j].V[k];
                statetemp.X[i][j] += W[i][j].V[k] * W[i][j].X[k];
                statetemp.Y[i][j] += W[i][j].V[k] * W[i][j].Y[k];
            }
        }
    }
    /*update intensities molecule by molecule*/
    double shrinkinv;
    double moltotalint;
    double moltotalintnew;
    double mollambda,moldeltaI;
    double x;
    for(i=0;i<statethis->N;i++){
	if(statethis->sp[i]=='1'){
	    shrinkinv = 1.0/(1.0+statethis->lambda[i]);
	    Intensity_Update_Homotopy(statethis->n[i], statetemp.I[i], run->alpha, statethis->I[i]);
	    for(j=0;j<statethis->n[i];j++){
		statethis->I[i][j] *= shrinkinv;
		if(statethis->I[i][j] < TrackingIntensityZero){
		    statethis->I[i][j] = TrackingIntensityZero;
		}
	    }
	}
	else if(statethis->sp[i]=='s'){
	    moltotalint = 0.0;
	    for(j=0;j<statethis->n[i];j++){
		moltotalint += statetemp.I[i][j];
	    }
	    mollambda = statethis->lambda[i];
	    moldeltaI = statethis->deltaI[i];
	    x = mollambda+moldeltaI-moltotalint;
            moltotalintnew = -x + sqrt(x*x+4.0*moldeltaI*moltotalint);
	    moltotalintnew = moltotalintnew*0.5;
	    shrinkinv = moltotalintnew / moltotalint;
	    Intensity_Update_Homotopy(statethis->n[i], statetemp.I[i], run->alpha, statethis->I[i]);
	    for(j=0;j<statethis->n[i];j++){
		statethis->I[i][j] *= shrinkinv;
		if(statethis->I[i][j] < TrackingIntensityZero){
		    statethis->I[i][j] = TrackingIntensityZero;
		}
	    }
	}
	else{
	}
    }
    /*update positions molecule by molecule*/
    double lambdainv;
    lambdainv = (run->sig * run->sig) / (run->D * run->D);
    for(i=0;i<statethis->N;i++){
	Inverse_EM_Matrix(statethis->n[i], statetemp.I[i], lambdainv, Inv_Matrix+i);
	Multiply_Matrix_Vector(Inv_Matrix+i, statetemp.X[i], statethis->X[i]);
	Multiply_Matrix_Vector(Inv_Matrix+i, statetemp.Y[i], statethis->Y[i]);
    }
/*compute log likelihood value*/
    double ans;
    /*imaging likelihood*/
    /*
    double sumintensity;
    ans = 0;
    sumintensity = 0;
    for(i=0;i<mv->Tmax;i++){
	for(j=0;j<sbs;j++){
	    if(mv->bound[i][j]){
	    }
	    else{
		ans += bb[i][j]*log(sumpsf[i][j]) - (1+run->lambda)*sumpsf[i][j] + (run->lambda)*Wn[i][j];
	    }
	}
    }
    */
    /*linking likelihood*/
    /*
    double Dinv;
    Dinv = 0.5/run->D;
    double logconst;
    logconst = log(2*M_PI*run->D*run->D);
    for(i=0;i<statethis->N;i++){
        for(j=0;j<statethis->n[i]-1;j++){
	    ans -= ((statethis->X[i][j]-statethis->X[i][j+1]) * (statethis->X[i][j]-statethis->X[i][j+1]) + (statethis->Y[i][j]-statethis->Y[i][j+1]) * (statethis->Y[i][j]-statethis->Y[i][j+1])) * Dinv;
            ans -= logconst;
        }
    }
    */
/*free malloced space*/
    for(i=0;i<mv->Tmax;i++){
	free(sumpsf[i]);
	free(sumpsfinv[i]);
	free(bb[i]);
    }
    free(bb);
    free(sumpsf);
    free(sumpsfinv);
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i];j++){
	    Sparse_free(W[i]+j);
	}
	free(W[i]);
    }
    free(W);
    for(i=0;i<mv->Tmax;i++){
	free(Wn[i]);
    }
    free(Wn);
    Free_State_Trace(&statetemp);
    for(i=0;i<statethis->N;i++){
        Free_EM_Matrix(Inv_Matrix+i);
    }
    free(Inv_Matrix);
/*return log likelihood*/
    return ans;
}


double EM_Tuning_LogLikelihood(struct EM_Trace_Run *run, struct Img_Series *mv, struct State_Trace *statethis, bool fix_frame0)
{
    double ans;
    double lossD=0.0;
    double lossI=0.0;
    double lossSP=0.0;
    double lossIMG=0.0;
    double lossInt=0.0;
    int i,j,k;
/*molecule motion loss*/
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i]-1;j++){
            lossD -= (statethis->X[i][j] - statethis->X[i][j+1])*(statethis->X[i][j] - statethis->X[i][j+1]);
            lossD -= (statethis->Y[i][j] - statethis->Y[i][j+1])*(statethis->Y[i][j] - statethis->Y[i][j+1]);
	}
    }
    lossD = 0.5*lossD/(run->D*run->D);
/*molecule intensity loss*/    
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i]-1;j++){
	    lossI -= run->alpha * fabs(log(statethis->I[i][j])-log(statethis->I[i][j+1]));
	}
    }
/*sparse penalty loss*/
    double molint;
    for(i=0;i<statethis->N;i++){
	molint = 0.0;
	for(j=0;j<statethis->n[i];j++){
	    molint += statethis->I[i][j];
	}
	if(statethis->sp[i]=='1'){
	    lossSP -= statethis->lambda[i]*molint;
	}
	else{
	    lossSP -= statethis->lambda[i]*log(1+molint/statethis->deltaI[i]);
	}
    }
/*Intensity loss*/
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i];j++){
            lossInt -= statethis->I[i][j];
	}
    }
/*observation loss*/
    int s1;
    int s2;
    s1 = run->s1;
    s2 = run->s2;
    int sb1;
    int sb2;
    sb1 = s1 + run->bsize + run->bsize;
    sb2 = s2 + run->bsize + run->bsize;
    int sbs;
    sbs = sb1*sb2;
    /*save the images, the boundary value of bb well be changed, we don't want to mess up the value inside run->bb */
    double **bb;
    bb = (double **)malloc(sizeof(double*)*(mv->Tmax));
    for(i=0;i<mv->Tmax;i++){
	bb[i] = (double*)malloc(sizeof(double)*sbs);
    }
    for(i=0;i<mv->Tmax;i++){
	for(j=0;j<sbs;j++){
	    bb[i][j] = mv->bb[i][j];
	}
    }
    /*save psf of the molecules*/
    struct Sparse **W;
    W = (struct Sparse **)malloc(sizeof(struct Sparse *)*(statethis->N));
    for(i=0;i<statethis->N;i++){
	W[i] = (struct Sparse *)malloc(sizeof(struct Sparse)*(statethis->n[i]));
    }
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i];j++){
	    Sparse_malloc(W[i]+j,(2+run->psfdecay+run->psfdecay)*(2+run->psfdecay+run->psfdecay));
	    W[i][j].n = 0;
	}
    }
    /*save psfno*/
    double **Wn;
    Wn = (double **)malloc(sizeof(double*)*(mv->Tmax));
    for(i=0;i<mv->Tmax;i++){
	Wn[i] = (double *)malloc(sizeof(double)*sbs);
    }
    /*save the sum up of psfs, expected images*/
    double **sumpsf;
    sumpsf = (double **)malloc(sizeof(double*)*(mv->Tmax));
    for(i=0;i<mv->Tmax;i++){
	sumpsf[i] = (double *)malloc(sizeof(double)*sbs);
    }
    /*compute and sum up psfs*/
    int idx;
    int T;
    /*compute psfs*/
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i];j++){
	    PSFGauss_Native(W[i]+j, run->exptable, run->expx, 0, sb1-1, 0, sb2-1, statethis->X[i][j], statethis->Y[i][j], run->sig, run->psfdecay);
	    for(k=0;k<W[i][j].n;k++){
	        W[i][j].V[k] = W[i][j].V[k] * statethis->I[i][j];
	    }
	}
    }
    /*compute psfno*/
    for(i=0;i<mv->Tmax;i++){
	PSFno(Wn[i],sbs,statethis->no[i]);
    }
    /*sum up psfs*/
    for(i=0;i<mv->Tmax;i++){
	for(j=0;j<sbs;j++){
	    sumpsf[i][j] = Wn[i][j];
	}
    }
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i];j++){
	    T = statethis->T[i][j];
	    for(k=0;k<W[i][j].n;k++){
		idx = W[i][j].X[k] + sb1*W[i][j].Y[k];
		sumpsf[T][idx] += W[i][j].V[k];
	    }
	}
    }
    /*compute the observation loss*/
    for(i=0;i<mv->Tmax;i++){
	for(j=0;j<sbs;j++){
	    if(mv->bound[i][j]){
	    }
	    else{
	        lossIMG += mv->bb[i][j]*log(sumpsf[i][j]);
	    }
	}
    }
/*free malloced space*/
    for(i=0;i<mv->Tmax;i++){
	free(sumpsf[i]);
	free(bb[i]);
    }
    free(bb);
    free(sumpsf);
    for(i=0;i<statethis->N;i++){
	for(j=0;j<statethis->n[i];j++){
	    Sparse_free(W[i]+j);
	}
	free(W[i]);
    }
    free(W);
    for(i=0;i<mv->Tmax;i++){
	free(Wn[i]);
    }
    free(Wn);

/*sum up losses and output*/
    ans = lossD+lossI+lossSP+lossIMG+lossInt;
    return ans;
}


