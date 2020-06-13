#ifndef ANALISH
#define ANALISH
extern double ystart[];
class Analis {
double xi, xi0, chimin, xmin, tmin;
const double fxi2;
void descent(double *py,double nv2[],int par[], int npf);
void revsvd(Mat_I_DP &aa,Vec_I_DP &w, Mat_I_DP &ai,Mat_O_DP &u);
void matmult(Mat_I_DP &aa, Mat_I_DP &ai,Mat_O_DP &u);
void hessian(int ndim,int np, Mat_DP& aa,Vec_DP& b,int par1[],double tmax);
public:
void sensitiv(const double tmax);
void coord(double ystart[]);
void grad(double tmax);
void sensitivity(Vec_DP &ystart);
Analis():fxi2(5.0){}
~Analis(){}
};
#endif

