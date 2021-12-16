#ifndef TrackingFunctions_H
#define TrackingFunctions_H
/*the index of a pixel in an image is x+y*s1*/
#include<math.h>
#include<stdlib.h>
#include"EMfunctions.h"
/*#include"MoleculeRandomize.h"*/

#define FactorialMax 5000
#define EM_matrix_zero 1e-14
#define EM_homotopy_zero 1e-14

/*v[0--n1-1][0--n2-1]*/
struct EM_Matrix
{
    int n1;
    int n2;
    double **v;
};

/*struct for intensity update, linked chain*/
struct Intensity_Segment
{
    int s; /*start of this segment*/
    int t; /*end of this segment*/
    char sign; /*larger or smaller or equal to next segment*/
    double I; /*Intensity*/
    struct Intensity_Segment *next;
};

/*molecule list of matrix form, [x,y,I,t,id]*/
struct State_Trace_Origin
{
    int l; /*l lines*/
    int Tmax; /*Tmax frames*/
    int IDnum;
    double *X;
    double *Y;
    double *I;
    int *T;
    int *ID;
    double *no;
};

/*Sophisticated representation of a tracking state*/
struct State_Trace
{
    int L; /*L frames*/
    int N; /*N molecule tracks*/
    int *l; /*l[i] is the number of molecules of frame i*/
    int *n; /*n[i] is the number of frames molecule i appears*/
    char *sp;  /*type of sparse penalty, '1' is L1, 's' is LSP, array sizeof of N*/
    double *lambda; /*sparse penalty, L1: lambda[i]*I, LSP: lambda[i]*log(1+I/deltaI[i]), array size of N*/
    double *deltaI; /*sparse penalty, LSP, array size of N*/
    int **idx1; /*X[idx1[t][0--l[t]-1]][idx2[t][0--l[t]-1]] enumerate all the X positions of the molecules in frame t*/
    int **idx2;
    double *no; /*noise intensity to L frames*/
    double **X; /*X[i][j]: i-->i'th molecule, j-->j'th time point since this molecule appeared */
    double **Y;
    double **I;
    int **T;
};

struct EM_Trace_Run
{
    int s1; /*lenth and width of the images*/
    int s2;
    int T; /*number of frames*/
    int psfdecay; /*longer than which psf decay to zero*/
    int bsize; /*extra boundary size*/
    int ITmax; /*can not exceed ITmax number of iterations*/
    double sig; /*width of molecule psf*/
    double D; /*expected diffusion velocity, Gaussian sigma*/
    double alpha; /*log Laplace parameter for intensity connection*/
/*in order not to compute exponential function, which is very time comsuming, we make it a table. exptable[i] = exp(expx[i])*/
    double *exptable; 
    double *expx;
};

struct SMC_Para
{
    /*number of particles*/
    int n;
    double B; /*Boltzmann factor*/
    double beta; /*energy = energyobserve+beta*energylinking*/
    double S; /*expected survival length(frames) of a molecule*/
    double A; /*expected activation rate, number per pixel area per frame*/
    double S0; /*for the first frame*/
    double A0; /*for the first frame*/
    double D; /*expected diffusion velocity, Gaussian sigma*/
    double alpha; /*log Laplace parameter for intensity connection*/
    double lambda; /*sparse penalty*/
    double *logfactable;
    double *exptable;
    double *expx;
};

struct SMC_Run
{
    int n;
    int t;
    struct State_Trace_Origin *particles;
    double *energy; /*log likelihood, probability proportional to exp(energy*BoltzmannFactor)*/
    double *eo; /*energyobserve, energy=energyobserve+para->beta*energylinking*/
    double *el; /*energylinking*/
};

struct Img_Series
{
    int s1;
    int s2;
    int Tmax;
    /*Extra boundary size*/
    int bsize;
    /*Image data, b[i] is the i'th frame*/
    double **b;
    /*Image data with extra boundary, pixels in extra boundary are assigned by zero*/
    double **bb;
    bool **bound; /*belong to boundary or not*/
};

struct EMTrack_Para{
    double sig;
    int bsize;
    int psfdecay;
    int maxite;
};

/*Usage example: QuickSortDouble(A,[0:n-1],0,n-1)*/
void QuickSortDouble(double *A, int *idx, int m, int n);
/*Usage example: QuickSortInt(A,[0:n-1],0,n-1)*/
void QuickSortInt(int *A, int *idx, int m, int n);

/*p->v[0--n1-1][0--n2-1]*/
void Malloc_EM_Matrix(struct EM_Matrix *p, int n1, int n2);
void Free_EM_Matrix(struct EM_Matrix *p);
/*Inverse of matrix for EM tracking, the matrix is of size n*n, the diagonal is "a+2*lambda(or a+lambda for the first or the last element)", v[i+1][i] = v[i][i+1] = -lambda*/
/*ans double be malloced(n*n) before calling*/
void Inverse_EM_Matrix(int n, double *a, double lambda, struct EM_Matrix *ans);

/*ans size p->n1, malloced before calling, v size p->n2*/
void Multiply_Matrix_Vector(struct EM_Matrix *p, double *v, double *ans);

/*free the space malloced inside pointer p*/
void Free_State_Trace_Origin(struct State_Trace_Origin *p);
/*re-assign the ID number in State_Trace_Origin to 0--IDnum-1, return IDnum*/
int Reassign_State_Trace_Origin(struct State_Trace_Origin *p);
/*Sort lines of State_Trace_Origin by time*/
void Sort_State_Trace_Origin(struct State_Trace_Origin *p);
/*malloc trace matrix and noise vector*/
void Malloc_State_Trace_Origin(struct State_Trace_Origin *p, int l, int Tmax);
/*po should be malloced before calling*/
void Deepcopy_State_Trace_Origin(struct State_Trace_Origin *pi, struct State_Trace_Origin *po);

/*convert struct State_Trace_Origin to struct State_Trace, re-assign of pi is required. pi->T should be sorted*/
void Make_State_Trace(struct State_Trace_Origin *pi, struct State_Trace *po);
/*malloc po as the same size to pi*/
void Malloc_State_Trace(struct State_Trace *pi, struct State_Trace *po);
/*po should be malloced before the calling*/
void Deepcopy_State_Trace(struct State_Trace *pi, struct State_Trace *po);
/*free the space malloced inside pointer p*/
void Free_State_Trace(struct State_Trace *p);

/*bb and bound are malloced inside*/
void Make_ExtraBoundary(struct Img_Series *p);

void Free_Img_Series(struct Img_Series *p);

/*ans should be malloced before calling, size n*/
void Intensity_Update_Homotopy(int n, double *a, double alpha, double *ans);



#endif


