#include"MoleculeRandomize.h"

void Mol_Survival_Randomize(int n, double svp, char *outans, std::mt19937_64 &gen)
{
    std::bernoulli_distribution bn_dis(svp);
    int i;
    for(i=0;i<n;i++){
	if(bn_dis(gen)){
	    outans[i] = 1;
	}
	else{
	    outans[i] = 0;
	}
    }
    return;
}


void Mol_Position_Randomize(int s1, int s2, double *b, int n, double *ansX, double *ansY, std::mt19937_64 &gen)
{
    int i;
    int ss=s1*s2;
    std::vector<double> fre(ss,0.0);
    for(i=0;i<ss;i++){
	fre[i] = b[i];
    }
    std::discrete_distribution<int> New_Mol_Gen(fre.begin(),fre.end());
    std::uniform_real_distribution<double> Mol_Pos_Refine_Gen(-0.50,0.50);
    int k;
    for(i=0;i<n;i++){
	k = New_Mol_Gen(gen);
	ansX[i] = k%s1+Mol_Pos_Refine_Gen(gen);
	ansY[i] = k/s1+Mol_Pos_Refine_Gen(gen);
    }
    return;
}


