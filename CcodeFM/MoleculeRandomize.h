#ifndef MoleculeRandomize_H
#define MoleculeRandomize_H
/*the index of a pixel in an image is x+y*s1*/
#include<stdio.h>
#include<stdlib.h>
#include<random>
#include<vector>
#include<chrono>

/*outans should be malloced outside of this function*/
void Mol_Survival_Randomize(int n, double svp, char *outans, std::mt19937_64 &gen);
/*ansX and ansY should be malloced outside of this function*/
/*b is assumed to be all non-negative*/
void Mol_Position_Randomize(int s1, int s2, double *b, int n, double *ansX, double *ansY, std::mt19937_64 &gen);
/*ansI should be malloced outside of this function*/
void Mol_Intensity_Randomize(int n, double *ansI, std::mt19937_64 &gen);

#endif
