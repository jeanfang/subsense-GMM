#pragma once
#include<algorithm>
#include<cmath>
#include<vector>
using namespace std;

class RTSuperpixel {
public:
	RTSuperpixel();
	virtual ~RTSuperpixel();
	void DBscan(unsigned char* R, unsigned char* G, unsigned char* B, int Rows, int Cols, unsigned short* label, double supnumber, int &realnumber, int post, unsigned int *labnumb);
	void mergin(int n, unsigned char* R, unsigned char* G, unsigned char* B, int Rows,
		int Cols, unsigned short* label, double *labelsL, double *labelsa, double *labelsb, double*labelsx, double *labelsy);
	int numb(int *a);
	void pixelQuery(int *neighbours, int neighboursP[5], int &num_neighb);
	void regionQuery(int* neighbours, int &k, unsigned char* R, unsigned char* G, unsigned  char* B, int n0, int n, int rows, int cols, unsigned short *label, bool *label0);
	void neighb0(unsigned short* label, int n, int Rows, int Cols, int ngb[5]);
	void supiel_neighbs(unsigned short* label, int *neighb, int ngb[5], int np, int Rows, int Cols, int &num_neighb);
};