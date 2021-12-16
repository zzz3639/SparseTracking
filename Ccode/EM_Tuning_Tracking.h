#ifndef EM_Tuning_Tracking_H
#define EM_Tuning_Tracking_H

#include<math.h>
#include<stdlib.h>
#include"EMfunctions.h"
#include"TrackingFunctions.h"

/*intensity of one molecule one frame can not below this value, to avoid divided by zero*/
#define TrackingIntensityZero 0.01

double EM_Tuning(struct Img_Series *mv, double sig, int psfdecay, int bsize, double D, double alpha, double lambda, struct State_Trace *state0, int ITmax, double llhzero, int llhwindow, bool fix_frame0, double *exptable, double *expx);

double EM_Tuning_Step(struct EM_Trace_Run *run, struct Img_Series *mv, struct State_Trace *statethis, bool fix_frame0);

double EM_Tuning_LogLikelihood(struct EM_Trace_Run *run, struct Img_Series *mv, struct State_Trace *statethis, bool fix_frame0);


#endif

