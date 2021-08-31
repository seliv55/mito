//---------------------------------------------------------------------------
#ifndef mainH
extern DP dxsav,*yy0, *yy1,xi0;  
extern Vec_INT *ija_p;
extern Vec_DP *xp_p,*sa_p,*x_p;
extern Mat_DP *yp_p,*yp_p,*d_p;
extern int nrhs;   // counts function evaluations
extern time_t ts,tf,tcal,tfirst;
extern void derivsl(const DP x, Vec_IO_DP &y, Vec_O_DP &dydx);
extern double integrbs(const double tmax,Vec_DP &ystart,const int);
double tisolve(const double tmax,double *y,int);
extern void pardep(Vec_DP &ystart,float ,bool flg=false);
extern double comb(double *py, int wr=0);

extern "C" void d02mwf_(int& NEQ,int& MAXORD,char& JCEVAL,double& HMAX,double& H0,int& ITOL,int ICOM[],int& LICOM,double COM[],int& LCOM,int& IFAIL);
extern "C" void d02nef_(int& NEQ,double& T,double& TOUT,double Y[],double YDOT[],double *RTOL,double *ATOL,int& ITASK,
void (*RES)(int& NEQ,double& T,double Y[],double YDOT[],double R[],int& IRES,int *IUSER,double *RUSER),
void (*JAC)(int NEQ,double T,double Y[],double YDOT[],double **PD,double& CJ,int *IUSER,double *RUSER),
int ICOM[],double COM[],int& LCOM,int *IUSER,double *RUSER,int& IFAIL);

extern "C" void resd02(int& NEQ,double& T,double Y[],double YDOT[],double R[],int& IRES,int *IUSER,double *RUSER);
extern void jacd02(int NEQ,double T,double Y[],double YDOT[],double **PD,double& CJ,int *IUSER,double *RUSER);
extern "C" void deliv_();

extern "C" void ddassl_(
void (*funcptr)(const double& time,  double y[], const double yprime[], double delta[], int& iRes, const double rPar[], const int iPar[]),
const int& noOfEquations,
double& currentTime,
double initialY[],
double initialYPrime[],
const double& finalTime,
int info[15],
double& relativeTolerance,
double& absoluteTolerance,
int& outputStatusFlag,
double dWorkArray[],
const int& lengthOfDWork,
int iWorkArray[],
const int& lengthOfIWork,
const double rParArray[],
const int iParArray[]
,void (*jacobian)(const double& time,  double y[], const double yprime[], double PD[][3], double& CJ,const double rPar[],const int iPar[])
);
    extern void res(const double& T, double y[],const double yprime[], double delta[], int& iRes, const double rPar[], const int iPar[]);
    extern void jac(const double& time,  double y[], const double yprime[], double dfdy[][3], double& CJ,const double rPar[],const  int iPar[]);
extern double d02nefsv(const double tout, double *py,int iter,int crv,double ros[]);
 
#define mainH
#endif
