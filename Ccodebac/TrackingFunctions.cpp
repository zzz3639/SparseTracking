#include"TrackingFunctions.h"
#include<vector>


void SwapDouble(double *x, double *y)
{
    double temp;
    temp = *x;
    *x = *y;
    *y = temp;
    return;
}

void SwapInt(int *x, int *y)
{
    int temp;
    temp = *x;
    *x = *y;
    *y = temp;
    return;
}


int ChoosePivot(int i,int j )
{
   return ((i+j)/2);
}

void QuickSortDouble(double *A, int *idx, int m, int n)
{
    double key;
    int i,j,k;
    if(m<n){
	k = ChoosePivot(m,n);
	SwapDouble(A+m,A+k);
	SwapInt(idx+m,idx+k);
	key = A[m];
	i = m+1;
	j = n;
	while(i<=j){
	    while((i<=n)&&(A[i]<=key))
		i++;
	    while((j>=m)&&(A[j]>key))
		j--;
	    if(i<j){
		SwapDouble(A+i,A+j);
		SwapInt(idx+i,idx+j);
	    }
	}
	SwapDouble(A+m,A+j);
	SwapInt(idx+m,idx+j);
	QuickSortDouble(A,idx,m,j-1);
	QuickSortDouble(A,idx,j+1,n);
    }
    return;
}

void QuickSortInt(int *A, int *idx, int m, int n)
{
    int key;
    int i,j,k;
    if(m<n){
	k = ChoosePivot(m,n);
	SwapInt(A+m,A+k);
	SwapInt(idx+m,idx+k);
	key = A[m];
	i = m+1;
	j = n;
	while(i<=j){
	    while((i<=n)&&(A[i]<=key))
		i++;
	    while((j>=m)&&(A[j]>key))
		j--;
	    if(i<j){
		SwapInt(A+i,A+j);
		SwapInt(idx+i,idx+j);
	    }
	}
	SwapInt(A+m,A+j);
	SwapInt(idx+m,idx+j);
	QuickSortInt(A,idx,m,j-1);
	QuickSortInt(A,idx,j+1,n);
    }
    return;
}


void Malloc_EM_Matrix(struct EM_Matrix *p, int n1, int n2)
{
    int i;
    p->n1 = n1;
    p->n2 = n2;
    p->v = (double**)malloc(sizeof(double*)*n1);
    for(i=0;i<n1;i++){
	p->v[i] = (double*)malloc(sizeof(double)*n2);
    }
    return;
}

void Free_EM_Matrix(struct EM_Matrix *p)
{
    int i;
    for(i=0;i<p->n1;i++){
	free(p->v[i]);
    }
    free(p->v);
    return;
}

void Inverse_EM_Matrix(int n, double *a, double lambda, struct EM_Matrix *ans)
{
    int i,j;
    if(n==0){
	return;
    }
    if(n==1){
	ans->v[0][0] = 1.0/a[0];
	return;
    }
    if(fabs(lambda)<EM_matrix_zero){
	for(i=0;i<n;i++){
	    ans->v[i][i] = 1.0/a[i];
	}
	return;
    }
    struct EM_Matrix M;
    Malloc_EM_Matrix(&M,n,n);
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            ans->v[i][j] = 0.0;
            M.v[i][j] = 0.0;
        }
	ans->v[i][i] = 1.0;
    }
    M.v[0][0] = a[0]+lambda;
    for(i=1;i<n-1;i++){
	M.v[i][i] = a[i]+2.0*lambda;
    }
    M.v[n-1][n-1] = a[n-1]+lambda;
    for(i=0;i<n-1;i++){
        M.v[i][i+1] = -lambda;
        M.v[i+1][i] = -lambda;
    }
    double x;
    for(i=0;i<n-1;i++){
        if(fabs(M.v[i][i])<EM_matrix_zero){
	    for(j=0;j<n;j++){
		M.v[i][j] += M.v[i+1][j];
		ans->v[i][j] += ans->v[i+1][j];
	    }
	}
	x = M.v[i+1][i]/M.v[i][i];
	for(j=0;j<n;j++){
	    M.v[i+1][j] -= x*M.v[i][j];
	    ans->v[i+1][j] -= x*ans->v[i][j];
	}
    }
    for(i=n-2;i>-1;i--){
	x = M.v[i][i+1]/M.v[i+1][i+1];
	for(j=0;j<n;j++){
	    M.v[i][j] -= x*M.v[i+1][j];
	    ans->v[i][j] -= x*ans->v[i+1][j];
	}
    }
    for(i=0;i<n;i++){
        x = 1.0/M.v[i][i];
	for(j=0;j<n;j++){
	    ans->v[i][j] *= x;
	}
    }
    Free_EM_Matrix(&M);
    return;
}

void Multiply_Matrix_Vector(struct EM_Matrix *p, double *v, double *ans)
{
    int i,j;
    for(i=0;i<p->n1;i++){
	ans[i] = 0;
    }
    for(i=0;i<p->n1;i++){
	for(j=0;j<p->n2;j++){
	    ans[i] += p->v[i][j]*v[j];
	}
    }
    return;
}


void Free_State_Trace_Origin(struct State_Trace_Origin *p)
{
    free(p->X);
    free(p->Y);
    free(p->I);
    free(p->T);
    free(p->ID);
    free(p->no);
    return;
}

void Free_State_Trace(struct State_Trace *p)
{
    int i;
    for(i=0;i<p->N;i++){
	free(p->X[i]);
	free(p->Y[i]);
	free(p->I[i]);
	free(p->T[i]);
    }
    free(p->X);
    free(p->Y);
    free(p->I);
    free(p->T);
    for(i=0;i<p->L;i++){
	free(p->idx1[i]);
	free(p->idx2[i]);
    }
    free(p->idx1);
    free(p->idx2);
    free(p->l);
    free(p->n);
    free(p->no);
    free(p->sp);
    free(p->lambda);
    free(p->deltaI);
    return;
}

void Sort_State_Trace_Origin(struct State_Trace_Origin *p)
{
    int i;
    int *intidx;
    intidx = (int*)malloc(sizeof(int)*p->l);
    for(i=0;i<p->l;i++){
	intidx[i] = i;
    }
    QuickSortInt(p->T,intidx,0,p->l-1);
    double *Xtemp;
    Xtemp = (double*)malloc(sizeof(double)*p->l);
    double *Ytemp;
    Ytemp = (double*)malloc(sizeof(double)*p->l);
    double *Itemp;
    Itemp = (double*)malloc(sizeof(double)*p->l);
    int *IDtemp;
    IDtemp = (int*)malloc(sizeof(int)*p->l);
    for(i=0;i<p->l;i++){
	Xtemp[i] = p->X[i];
	Ytemp[i] = p->Y[i];
	Itemp[i] = p->I[i];
	IDtemp[i] = p->ID[i];
    }
    for(i=0;i<p->l;i++){
	p->X[i] = Xtemp[intidx[i]];
	p->Y[i] = Ytemp[intidx[i]];
	p->I[i] = Itemp[intidx[i]];
	p->ID[i] = IDtemp[intidx[i]];
    }
    free(intidx);
    free(Xtemp);
    free(Ytemp);
    free(Itemp);
    free(IDtemp);
    return;
}

int Reassign_State_Trace_Origin(struct State_Trace_Origin *p)
{
    int IDmax=-1;
    int i;
    int *ID=p->ID;
    int l=p->l;
    for(i=0;i<l;i++){
	if(ID[i]>IDmax){
	    IDmax = ID[i];
	}
    }
    IDmax = IDmax+1;
    char *IDbarrel;
    IDbarrel = (char *)malloc(sizeof(char)*IDmax);
    for(i=0;i<IDmax;i++){
	IDbarrel[i] = 0;
    }
    for(i=0;i<l;i++){
	IDbarrel[ID[i]] = 1;
    }
    int *IDidx;
    int k=0;
    IDidx = (int*)malloc(sizeof(int)*IDmax);
    for(i=0;i<IDmax;i++){
	if(IDbarrel[i]==0){
	    continue;
	}
	IDidx[i] = k;
	k++;
    }
    p->IDnum = k;
    for(i=0;i<l;i++){
	ID[i] = IDidx[ID[i]];
    }
    free(IDbarrel);
    free(IDidx);
    return k;
}

void Make_State_Trace(struct State_Trace_Origin *pi, struct State_Trace *po)
{
    int IDnum;
    IDnum = pi->IDnum;
    int l=pi->l;
    int Tmax=pi->Tmax;
    po->L = Tmax;
    po->l = (int*)malloc(sizeof(int)*Tmax);
    po->N = IDnum;
    po->n = (int*)malloc(sizeof(int)*IDnum);
    po->no = (double*)malloc(sizeof(double)*Tmax);

    po->sp = (char*)malloc(sizeof(char)*IDnum);
    po->lambda = (double*)malloc(sizeof(double)*IDnum);
    po->deltaI = (double*)malloc(sizeof(double)*IDnum);

    int i,j;
    for(i=0;i<Tmax;i++){
        po->no[i] = pi->no[i];
	po->l[i] = 0;
    }
    for(i=0;i<IDnum;i++){
	po->n[i] = 0;
    }
    for(i=0;i<l;i++){
	po->l[pi->T[i]]++;
	po->n[pi->ID[i]]++;
    }
    po->idx1 = (int**)malloc(sizeof(int*)*Tmax);
    po->idx2 = (int**)malloc(sizeof(int*)*Tmax);
    for(i=0;i<Tmax;i++){
	po->idx1[i] = (int*)malloc(sizeof(int)*po->l[i]);
	po->idx2[i] = (int*)malloc(sizeof(int)*po->l[i]);
    }
    po->X = (double**)malloc(sizeof(double*)*IDnum);
    po->Y = (double**)malloc(sizeof(double*)*IDnum);
    po->I = (double**)malloc(sizeof(double*)*IDnum);
    po->T = (int**)malloc(sizeof(int*)*IDnum);
    for(i=0;i<IDnum;i++){
	po->X[i] = (double*)malloc(sizeof(double)*po->n[i]);
	po->Y[i] = (double*)malloc(sizeof(double)*po->n[i]);
	po->I[i] = (double*)malloc(sizeof(double)*po->n[i]);
	po->T[i] = (int*)malloc(sizeof(int)*po->n[i]);
    }
    int *IDarray;
    IDarray = (int*)malloc(sizeof(int)*IDnum);
    for(i=0;i<IDnum;i++){
	IDarray[i] = 0;
    }
    int *Tarray;
    Tarray = (int*)malloc(sizeof(int)*Tmax);
    for(i=0;i<Tmax;i++){
	Tarray[i] = 0;
    }
    for(i=0;i<l;i++){
        po->X[pi->ID[i]][IDarray[pi->ID[i]]] = pi->X[i];
        po->Y[pi->ID[i]][IDarray[pi->ID[i]]] = pi->Y[i];
        po->I[pi->ID[i]][IDarray[pi->ID[i]]] = pi->I[i];
        po->T[pi->ID[i]][IDarray[pi->ID[i]]] = pi->T[i];
	po->idx1[pi->T[i]][Tarray[pi->T[i]]] = pi->ID[i];
	po->idx2[pi->T[i]][Tarray[pi->T[i]]] = IDarray[pi->ID[i]];
	IDarray[pi->ID[i]]++;
	Tarray[pi->T[i]]++;
    }
    free(IDarray);
    free(Tarray);
    return;
}


void Malloc_State_Trace_Origin(struct State_Trace_Origin *p, int l, int Tmax)
{
    p->l = l;
    p->Tmax = Tmax;
    p->X = (double*)malloc(sizeof(double)*l);
    p->Y = (double*)malloc(sizeof(double)*l);
    p->I = (double*)malloc(sizeof(double)*l);
    p->T = (int*)malloc(sizeof(int)*l);
    p->ID= (int*)malloc(sizeof(int)*l);
    p->no= (double*)malloc(sizeof(double)*Tmax);
    return;
}

void Deepcopy_State_Trace_Origin(struct State_Trace_Origin *pi, struct State_Trace_Origin *po)
{
    int i;
    po->l = pi->l;
    po->Tmax = pi->Tmax;
    po->IDnum = pi->IDnum;
    for(i=0;i<pi->Tmax;i++){
	po->no[i] = pi->no[i];
    }
    for(i=0;i<pi->l;i++){
	po->X[i] = pi->X[i];
	po->Y[i] = pi->Y[i];
	po->I[i] = pi->I[i];
	po->T[i] = pi->T[i];
	po->ID[i]=pi->ID[i];
    }
    return;
}

void Malloc_State_Trace(struct State_Trace *pi, struct State_Trace *po)
{
    int i;
    po->L = pi->L;
    po->N = pi->N;
    po->l = (int*)malloc(sizeof(int)*(pi->L));
    po->n = (int*)malloc(sizeof(int)*(pi->N));
    po->no = (double*)malloc(sizeof(double)*(pi->L));

    po->sp = (char*)malloc(sizeof(char)*(pi->N));
    po->lambda = (double*)malloc(sizeof(double)*(pi->N));
    po->deltaI = (double*)malloc(sizeof(double)*(pi->N));

    po->idx1 = (int**)malloc(sizeof(int*)*pi->L);
    po->idx2 = (int**)malloc(sizeof(int*)*pi->L);
    for(i=0;i<pi->L;i++){
	po->idx1[i] = (int*)malloc(sizeof(int)*pi->l[i]);
	po->idx2[i] = (int*)malloc(sizeof(int)*pi->l[i]);
    }

    po->X = (double**)malloc(sizeof(double*)*pi->N);
    po->Y = (double**)malloc(sizeof(double*)*pi->N);
    po->I = (double**)malloc(sizeof(double*)*pi->N);
    po->T = (int**)malloc(sizeof(int*)*pi->N);
    for(i=0;i<pi->N;i++){
	po->X[i] = (double*)malloc(sizeof(double)*pi->n[i]);
	po->Y[i] = (double*)malloc(sizeof(double)*pi->n[i]);
        po->I[i] = (double*)malloc(sizeof(double)*pi->n[i]);
	po->T[i] = (int*)malloc(sizeof(int)*pi->n[i]);
    }
    return;
}

void Deepcopy_State_Trace(struct State_Trace *pi, struct State_Trace *po)
{
    int i,j;
    po->L = pi->L;
    po->N = pi->N;
    for(i=0;i<pi->L;i++){
	po->no[i] = pi->no[i];
	po->l[i] = pi->l[i];
	for(j=0;j<pi->l[i];j++){
	    po->idx1[i][j] = pi->idx1[i][j];
	    po->idx2[i][j] = pi->idx2[i][j];
	}
    }
    for(i=0;i<pi->N;i++){
	po->n[i] = pi->n[i];
        po->sp[i] = pi->sp[i];
        po->lambda[i] = pi->lambda[i];
        po->deltaI[i] = pi->deltaI[i];
	for(j=0;j<pi->n[i];j++){
	    po->X[i][j] = pi->X[i][j];
	    po->Y[i][j] = pi->Y[i][j];
	    po->I[i][j] = pi->I[i][j];
	    po->T[i][j] = pi->T[i][j];
	}
    }
    return;
}

void Make_ExtraBoundary(struct Img_Series *p)
{
    int i,j;
    int sb1=p->s1+p->bsize+p->bsize; 
    int sb2=p->s2+p->bsize+p->bsize;
    int sbs=sb1*sb2;
    p->bb = (double**)malloc(sizeof(double*)*p->Tmax);
    p->bound = (bool**)malloc(sizeof(bool*)*p->Tmax);
    for(i=0;i<p->Tmax;i++){
	p->bb[i] = (double*)malloc(sizeof(double)*sbs);
	p->bound[i] = (bool*)malloc(sizeof(bool)*sbs);
    }
    for(i=0;i<p->Tmax;i++){
	for(j=0;j<sbs;j++){
	    p->bb[i][j] = 0;
	    p->bound[i][j] = true;
	}
    }
    int a,b;
    int aa;
    for(i=0;i<p->Tmax;i++){
	for(a=0;a<p->s1;a++){
	    aa = a+p->bsize;
	    for(b=0;b<p->s2;b++){
		p->bb[i][aa+(b+p->bsize)*sb1] = p->b[i][a+b*p->s1];
		p->bound[i][aa+(b+p->bsize)*sb1] = false;
	    }
	}
    }
    return;
}

void Free_Img_Series(struct Img_Series *p)
{
    int i;
    for(i=0;i<p->Tmax;i++){
	free(p->b[i]);
	free(p->bb[i]);
	free(p->bound[i]);
    }
    free(p->b);
    free(p->bb);
    free(p->bound);
    return;
}

void Intensity_Update_Homotopy(int n, double *a, double alpha, double *ans)
{
/*initialize*/
    if(n<=0){
	return;
    }
    if(n==1){
	ans[0] = a[0];
	return;
    }
    int i;
    if(alpha<EM_homotopy_zero){
	for(i=0;i<n;i++){
	    ans[i] = a[i];
	}
	return;
    }
    double *epsilon;
    epsilon = (double*)malloc(sizeof(double)*(n-1));
    double x,y,z;
    int j,k;
    x = 0.0;
    for(i=0;i<n;i++){
	x += a[i];
    }
    x = x/n;
    epsilon[0] = a[0] - x;
    for(i=1;i<n-1;i++){
	epsilon[i] = epsilon[i-1] + a[i] - x;
    }
    double alphathis;
    alphathis = fabs(epsilon[0]);
    for(i=1;i<n-1;i++){
	y = fabs(epsilon[i]);
	if(y>alphathis){
	    alphathis = y;
	}
    }
    if(alpha>=alphathis){
	for(i=0;i<n;i++){
	    ans[i] = x;
	}
        free(epsilon);
	return;
    }
    alphathis = alphathis+1.0;
    struct Intensity_Segment segchain;
    segchain.next = (struct Intensity_Segment *)malloc(sizeof(struct Intensity_Segment));
    segchain.next->s = 0;
    segchain.next->t = n;
    segchain.next->sign = 0;
    segchain.next->I = x;
    segchain.next->next = NULL;
    struct Intensity_Segment *sp1;
    struct Intensity_Segment *sp2;
    struct Intensity_Segment *spnew;
    struct Intensity_Segment *spmax;
    double *depsilon;
    depsilon = (double*)malloc(sizeof(double)*(n-1));
/*loop of homotopy algorithm*/
    char signbefore,signafter;
    double alphabefore,alphaafter;
    double alphachange;
    double epsilonthis;
    double dalphathis;
    double alphanexttemp;
    double alphanexttempplus;
    double alphanexttempminus;
    double alphanextmax;
    int bpoint;
    char bpointc;
    while(1){
    /*judge if other break point exists*/
        sp1 = &segchain;
	sp2 = sp1->next;
	signbefore = 0;
	signafter = sp2->sign;
	for(i=0;i<n-1;i++){
	    epsilon[i] = 0.0;
	}
	bpointc = 0;
	while(sp2!=NULL){
	    if(sp2->s==sp2->t-1){
	    }
	    else{
                alphabefore = signbefore*alphathis;
                alphaafter = -signafter*alphathis;
		alphachange = alphabefore + alphaafter;
		x = alphachange;
		for(i=sp2->s;i<sp2->t;i++){
		    x += a[i];
		}
		x = x/(sp2->t-sp2->s);
		y = 1.0/(sp2->t-sp2->s);
		epsilon[sp2->s] = a[sp2->s] + alphabefore - x;
		depsilon[sp2->s]= signbefore - y*signbefore + y*signafter;
		for(i=sp2->s+1;i<sp2->t-2;i++){
		    epsilon[i] = epsilon[i-1] + a[i] - x;
		    depsilon[i] = depsilon[i-1] - y*signbefore + y*signafter;
		}
		epsilon[sp2->t-2] = x - a[sp2->t-1] - alphaafter;
		depsilon[sp2->t-2] = signafter + y*signbefore - y*signafter;
		bpointc = 0;
                for(i=sp2->s;i<sp2->t-1;i++){
		    if(alphathis-fabs(epsilon[i])<EM_homotopy_zero){
			bpoint = i;
			bpointc= 1;
			break;
		    }
		}
		if(bpointc){
		    break;
		}
	    }
	    sp1 = sp1->next;
	    sp2 = sp2->next;
	    signbefore = sp1->sign;
	    if(sp2==NULL){
		signafter = 0;
	    }
	    else{
		signafter = sp2->sign;
	    }
	}
    /*update*/
        if(bpointc){
	    spnew = (struct Intensity_Segment *)malloc(sizeof(struct Intensity_Segment));
	    spnew->s = bpoint + 1;
	    spnew->t = sp2->t;
	    spnew->sign = sp2->sign;
	    spnew->next = sp2->next;
	    sp2->t = bpoint + 1;
	    sp2->next = spnew;
	    if(epsilon[bpoint]>0){
		sp2->sign = 1;
	    }
	    else{
		sp2->sign = -1;
	    }
	    continue;
        }
    /*find next alpha and next break point*/
	alphanextmax = 0.0;
        sp1 = &segchain;
	sp2 = sp1->next;
	signbefore = 0;
	signafter = sp2->sign;
	for(i=0;i<n-1;i++){
	    epsilon[i] = 0.0;
	    depsilon[i] = 0.0;
	}
	while(sp2!=NULL){
	    if(sp2->s==sp2->t-1){
	    }
	    else{
                alphabefore = signbefore*alphathis;
                alphaafter = -signafter*alphathis;
		alphachange = alphabefore + alphaafter;
		x = alphachange;
		for(i=sp2->s;i<sp2->t;i++){
		    x += a[i];
		}
		x = x/(sp2->t-sp2->s);
		y = 1.0/(sp2->t-sp2->s);
		epsilon[sp2->s] = a[sp2->s] + alphabefore - x;
		depsilon[sp2->s]= signbefore - y*signbefore + y*signafter;
		for(i=sp2->s+1;i<sp2->t-2;i++){
		    epsilon[i] = epsilon[i-1] + a[i] - x;
		    depsilon[i] = depsilon[i-1] - y*signbefore + y*signafter;
		}
		epsilon[sp2->t-2] = x - a[sp2->t-1] - alphaafter;
		depsilon[sp2->t-2] = signafter + y*signbefore - y*signafter;
		for(i=sp2->s;i<sp2->t-1;i++){
                /*count alphanexttemp*/
		    x = alphathis - epsilon[i];
		    y = 1.0 - depsilon[i];
                    if(y<EM_homotopy_zero){
			alphanexttempplus = 0.0;
		    }
		    else{
			z = alphathis -  x/y;
			alphanexttempplus = MAX(0.0,z);
		    }
		    x = alphathis + epsilon[i];
		    y = 1.0 + depsilon[i];
		    if(y<EM_homotopy_zero){
			alphanexttempminus = 0.0;
		    }
		    else{
			z = alphathis - x/y;
			alphanexttempminus = MAX(0.0,z);
		    }
		    alphanexttemp = MAX(alphanexttempplus,alphanexttempminus);
		    if(alphanexttemp>alphanextmax){
			alphanextmax = alphanexttemp;
			bpoint = i;
			spmax = sp2;
		    }
		}
	    }
	    sp1 = sp1->next;
	    sp2 = sp2->next;
	    signbefore = sp1->sign;
	    if(sp2==NULL){
		signafter = 0;
	    }
	    else{
		signafter = sp2->sign;
	    }
	}
    /*update*/
	if(alphanextmax<EM_homotopy_zero){
	    break;
	}
	if(alphanextmax<alpha){
	    break;
	}
	spnew = (struct Intensity_Segment *)malloc(sizeof(struct Intensity_Segment));
	spnew->s = bpoint + 1;
	spnew->t = spmax->t;
	spnew->sign = spmax->sign;
	spnew->next = spmax->next;
	spmax->t = bpoint + 1;
	spmax->next = spnew;
	if(epsilon[bpoint]>0){
	    spmax->sign = 1;
	}
	else{
	    spmax->sign = -1;
	}
	alphathis = alphanextmax;
    }
/*count final result and output to ans*/

    sp2 = segchain.next;
    signbefore = 0;
    signafter = sp2->sign;
    while(sp2!=NULL){
	alphabefore = signbefore*alpha;
	alphaafter = -signafter*alpha;
	alphachange = alphabefore + alphaafter;
	x = alphachange;
	for(i=sp2->s;i<sp2->t;i++){
	    x += a[i];
	}
	x = x/(sp2->t-sp2->s);
	for(i=sp2->s;i<sp2->t;i++){
	    ans[i] = x;
	}
	signbefore = sp2->sign;
	sp2 = sp2->next;
	if(sp2==NULL){
	    signafter = 0;
	}
	else{
	    signafter = sp2->sign;
	}
    }

/*free malloced space*/
    sp1 = segchain.next;
    while(sp1!=NULL){
	sp2 = sp1->next;
	free(sp1);
	sp1 = sp2;
    }
    free(depsilon);
    free(epsilon);
    return;
}

void SMC_Initialize(struct SMC_Run *r, struct State *pic0, struct Img_Series *img_full, struct SMC_Para *para, struct EMImg_Para *empara, std::mt19937_64 &gen)
{
    int ss;
    ss = img_full->s1*img_full->s2;
    double *resimg;
    resimg = (double*)malloc(sizeof(double)*ss);
    int i;
    double meanI = 0.0;
    for(i=0;i<pic0->n;i++){
	meanI += pic0->I[i];
    }
    meanI = meanI/pic0->n;
    r->n = para->n;
    r->t = 1;
    r->particles = (struct State_Trace_Origin*)malloc(sizeof(struct State_Trace_Origin)*para->n);
    r->energy = (double*)malloc(sizeof(double)*para->n);
    r->el = (double*)malloc(sizeof(double)*para->n);
    r->eo = (double*)malloc(sizeof(double)*para->n);
    struct State picthis;
    struct State pictemp;
    struct RUN runthis;
    int sbs;
    sbs = (img_full->s1+empara->psfdecay+empara->psfdecay)*(img_full->s2+empara->psfdecay+empara->psfdecay);
    runthis.bb = (double*)malloc(sizeof(double)*sbs);
    runthis.bound = (bool*)malloc(sizeof(bool)*sbs);
    runthis.Wn= (double*)malloc(sizeof(double)*sbs);
    runthis.S = (double*)malloc(sizeof(double)*sbs);
    for(i=0;i<sbs;i++){
	runthis.bb[i] = img_full->bb[0][i];
	runthis.bound[i] = img_full->bound[0][i];
    }
    int numac;
    int numsv;
    int nummol;
    char *svvec;
    svvec = (char*)malloc(sizeof(char)*r->n);
    double *ansX;
    double *ansY;
    double *ansI;
    double avac=para->A0*img_full->s1*img_full->s2;
    double svp=1.0-1.0/(para->S0+1.0);
    std::poisson_distribution<int> disac(avac);
    int j,k;
    double *expimgtemp;
    expimgtemp = (double*)malloc(sizeof(double)*ss);
    double eo,el;
    for(i=0;i<para->n;i++){
    /*do sampling*/
        numac = disac(gen);
	Mol_Survival_Randomize(pic0->n,svp,svvec,gen);
	numsv = 0;
	for(j=0;j<pic0->n;j++){
            if(svvec[j]){
		numsv++;
	    }
	}
        nummol = numac+numsv;
	ansX = (double*)malloc(sizeof(double)*numac);
	ansY = (double*)malloc(sizeof(double)*numac);
	ansI = (double*)malloc(sizeof(double)*numac);
	picthis.no = pic0->no;
	picthis.X = (double*)malloc(sizeof(double)*nummol);
	picthis.Y = (double*)malloc(sizeof(double)*nummol);
	picthis.I = (double*)malloc(sizeof(double)*nummol);
	k = 0;
	for(j=0;j<pic0->n;j++){
	    if(svvec[j]){
		picthis.X[k] = pic0->X[j];
		picthis.Y[k] = pic0->Y[j];
		picthis.I[k] = pic0->I[j];
		k++;
	    }
	}
	picthis.n = k;
	ExpectedImg(&picthis, empara->sig, img_full->s1, img_full->s2, empara->psfdecay, resimg, para->exptable, para->expx);
	for(j=0;j<ss;j++){
	    resimg[j] = fabs(resimg[j]-img_full->b[0][j]);
	}
	Mol_Position_Randomize(img_full->s1,img_full->s2,resimg,numac,ansX,ansY,gen);
        for(j=0;j<numac;j++){
	    ansI[j] = meanI;
	}
	for(j=0;j<numac;j++){
	    picthis.X[k+j] = ansX[j];
	    picthis.Y[k+j] = ansY[j];
	    picthis.I[k+j] = ansI[j];
	}
	picthis.n = nummol;
    /*tunning by EM algorithm*/
        /*Initialize struct RUN for EM algorithm*/
	State_malloc(&pictemp,picthis.n,0,true);
	runthis.n = picthis.n;
	runthis.sig = empara->sig;
	runthis.bsize = empara->bsize;
	runthis.psfdecay = empara->psfdecay;
	runthis.lambda = empara->lambda;
	runthis.Pzero = empara->Pzero;
	runthis.Izero = empara->Izero;
	runthis.s1 = img_full->s1;
	runthis.s2 = img_full->s2;
	runthis.sb1 = img_full->s1+empara->psfdecay+empara->psfdecay;
	runthis.sb2 = img_full->s2+empara->psfdecay+empara->psfdecay;
	runthis.ss = img_full->s1*img_full->s2;
	runthis.sbs = runthis.sb1*runthis.sb2;
	runthis.W = (struct Sparse*)malloc(sizeof(struct Sparse)*picthis.n);
	runthis.w = (struct Sparse*)malloc(sizeof(struct Sparse)*picthis.n);
	for(j=0;j<picthis.n;j++){
	    Sparse_malloc(runthis.W+j,(2+empara->psfdecay+empara->psfdecay)*(2+empara->psfdecay+empara->psfdecay));
	    Sparse_malloc(runthis.w+j,(2+empara->psfdecay+empara->psfdecay)*(2+empara->psfdecay+empara->psfdecay));
	    runthis.W[j].n=0;
	    runthis.w[j].n=0;
	}
	runthis.terminate = 'N';
        runthis.exptable = para->exptable;
	runthis.expx = para->expx;
	runthis.X = (double*)malloc(sizeof(double)*picthis.n);
	runthis.Y = (double*)malloc(sizeof(double)*picthis.n);
	runthis.die = (bool*)malloc(sizeof(bool)*picthis.n);
	runthis.freeze = (bool*)malloc(sizeof(bool)*picthis.n);
	for(j=0;j<picthis.n;j++){
	    runthis.X[j] = -1.0;
	    runthis.Y[j] = -1.0;
	    runthis.die[j] = false;
	    runthis.freeze[j] = false;
	}
        /*Run EM algorithm step by step until converge*/
	for(j=0;j<empara->maxite;j++){
	    RunStep(&runthis,&picthis,&pictemp,false,false,NULL,true);
	    if(runthis.terminate=='Y')
		break;
	}
    /*record result*/
	/*record molecule table*/
	Malloc_State_Trace_Origin(r->particles+i,picthis.n,1);
	for(j=0;j<picthis.n;j++){
	    r->particles[i].X[j] = picthis.X[j];
	    r->particles[i].Y[j] = picthis.Y[j];
	    r->particles[i].I[j] = picthis.I[j];
	    r->particles[i].T[j] = 0;
	    r->particles[i].ID[j] = j;
	    r->particles[i].IDnum = picthis.n;
	    r->particles[i].no[0] = picthis.no;
	}
	/*record energy*/
        ExpectedImg(&picthis, empara->sig, img_full->s1, img_full->s2, empara->psfdecay, expimgtemp, para->exptable, para->expx);
	eo = 0.0;
	for(j=0;j<ss;j++){
            eo -= img_full->b[0][j]*log(expimgtemp[j]) - expimgtemp[j];
	}
	el = 0.0;
	r->el[i] = el;
	r->eo[i] = eo;
	r->energy[i] = r->eo[i];
    /*free spaces*/
	State_free(&pictemp,true);
	free(runthis.X);
	free(runthis.Y);
	free(runthis.die);
	free(runthis.freeze);
	for(j=0;j<picthis.n;j++){
	    Sparse_free(runthis.W+j);
	    Sparse_free(runthis.w+j);
	}
	free(runthis.W);
	free(runthis.w);
	free(ansX);
	free(ansY);
	free(ansI);
	free(picthis.X);
	free(picthis.Y);
	free(picthis.I);
    }
    /*free spaces*/
    free(runthis.bb);
    free(runthis.bound);
    free(runthis.S);
    free(runthis.Wn);
    free(resimg);
    free(expimgtemp);
    return;
}

void SMC_Para_Initialize(struct SMC_Para *para)
{
    int i;
    para->exptable=(double *)malloc(sizeof(double)*Expmin*Expinterval);
    para->expx=(double *)malloc(sizeof(double)*Expmin*Expinterval);
    int k;
    double dx;
    k=Expmin*Expinterval;
    dx=1.0/Expinterval;
    for(i=0;i<k;i++){
	para->expx[i]=-dx*i;
	para->exptable[i]=exp(-dx*i);
    }
    para->logfactable=(double *)malloc(sizeof(double)*FactorialMax);
    para->logfactable[0] = 0;
    for(i=1;i<FactorialMax;i++){
	para->logfactable[i] = para->logfactable[i-1] + log((double)i);
    }
    return;
}

/*sequencial monte carlo, also known as particle filter*/
void SMC_Step(struct SMC_Run *r, struct Img_Series *img_full, struct SMC_Para *para, struct EMTrack_Para *emtpara, std::mt19937_64 &gen)
{
    double *ansX;
    double *ansY;
    double avac=para->A*img_full->s1*img_full->s2;
    double svp=1.0-1.0/(para->S+1.0);
    std::vector<double> probparticles(r->n,0.0);
    int i;
    for(i=0;i<r->n;i++){
	probparticles[i] = exp(r->energy[i]*para->B);
    }
    std::discrete_distribution<int> pselector(probparticles.begin(),probparticles.end());
    int h,j,k;
    int indp;
    struct State_Trace_Origin *topointer;
    struct State_Trace_Origin traceoriginnext;
    struct State_Trace tracenext;
    struct State picnext;
    int numsv,numac,nummol;
    std::poisson_distribution<int> disac(avac);
    char *svvec;
    double *resimg;
    int ss;
    ss = img_full->s1*img_full->s2;
    resimg = (double*)malloc(sizeof(double)*ss);
    double meanI;
    for(i=0;i<r->n;i++){
	indp = pselector(gen);
        topointer = r->particles+indp;
	k = Reassign_State_Trace_Origin(topointer);
	topointer->IDnum = k;
	/*propose new particles*/
	k=0;
        for(j=0;j<topointer->l;j++){
	    if(topointer->T[j]==r->t-1){
		k++;
	    }
	}
	svvec = (char*)malloc(sizeof(char)*k);
	Mol_Survival_Randomize(k,svp,svvec,gen);
	numsv = 0;
	for(j=0;j<k;j++){
	    if(svvec[j]){
		numsv++;
	    }
	}
        numac = disac(gen);
	nummol = numac+numsv;
	Malloc_State_Trace_Origin(&traceoriginnext, topointer->l+nummol, r->t+1);
	State_malloc(&picnext,numsv,0,true);
	Deepcopy_State_Trace_Origin(topointer,&traceoriginnext);
        k=0;
	h=0;
	picnext.n = numsv;
	picnext.no= topointer->no[r->t-1];
	meanI = 0.0;
	for(j=0;j<topointer->l;j++){
	    if(topointer->T[j]==r->t-1){
		meanI += topointer->I[j];
		if(svvec[k]){
		    picnext.X[h] = topointer->X[j];
		    picnext.Y[h] = topointer->Y[j];
                    picnext.I[h] = topointer->I[j];
		    traceoriginnext.X[topointer->l+h] = topointer->X[j];
		    traceoriginnext.Y[topointer->l+h] = topointer->Y[j];
		    traceoriginnext.I[topointer->l+h] = topointer->I[j];
		    traceoriginnext.T[topointer->l+h] = r->t;
		    traceoriginnext.ID[topointer->l+h]= topointer->ID[j];
		    h++;
		}
		k++;
	    }
	}
        meanI = meanI/k;
	ExpectedImg(&picnext, emtpara->sig, img_full->s1, img_full->s2, emtpara->psfdecay, resimg, para->exptable, para->expx);
	for(j=0;j<ss;j++){
	    resimg[j] = fabs(resimg[j]-img_full->b[r->t][j]);
	}
	ansX = (double*)malloc(sizeof(double)*numac);
	ansY = (double*)malloc(sizeof(double)*numac);
	Mol_Position_Randomize(img_full->s1,img_full->s2,resimg,numac,ansX,ansY,gen);
	for(j=0;j<numac;j++){
	    traceoriginnext.X[topointer->l+numsv+j] = ansX[j];
	    traceoriginnext.Y[topointer->l+numsv+j] = ansY[j];
	    traceoriginnext.I[topointer->l+numsv+j] = meanI;
	    traceoriginnext.T[topointer->l+numsv+j] = r->t;
	    traceoriginnext.ID[topointer->l+numsv+j]= topointer->IDnum+j;
	}
	traceoriginnext.IDnum = topointer->IDnum + nummol;
	/*Tuning the particles*/
	/*Update sampling probabilities(energy) of the particle set*/
	/*free spaces*/
	free(ansX);
	free(ansY);
	free(svvec);
	State_free(&picnext,true);
	Free_State_Trace_Origin(&traceoriginnext);
    }
    free(resimg);
    return;
}





