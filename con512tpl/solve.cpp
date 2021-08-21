#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm> 
#include "nr.h"
#include "StiffIntegratorT.h"
#include "NonStiffIntegratorT.h"
#include "nums.hh"
#include "modlab.h"
#include "main.hh"
//---------------------------------------------------------------------------
using namespace std;
using namespace Label;
  vector<double> vx;

void isores(const double& T, double y[],const double yprime[], double delta[], int& iRes, const double rPar[], const int iPar[]){
       double f[NN];
	horse.distr(y, f);
	for(unsigned i=0;i<NN;i++) delta[i]=f[i]-yprime[i];
}
void jac(const double& time,  double y[], const double yPrime[], double dfdy[][3], double& CJ, const double rPar[], const int iPar[]){}

int Ldistr::ddisolve(int istart,const double tfin, double *py,ostringstream& fkin) {
  int info[15],idid=0,lrw=200000,liw=500,iwork[500],ipar[2];
  double tinc=0.35*tfac,to,rtol=0.001,atol=1.0e-9,h0=1e-5, hmax=2.;
  double t=istart*tinc,term=tfin*tfac,rpar[2], yprime[NN], rwork[200000];
  rwork[1]=hmax; rwork[2]=h0;
  for(int i=0;i<15;i++) info[i]=0;   info[6]=1; info[10]=1;
  distr(py, yprime);
//   for(int i=0;i<NN;i++) cout<<i<<") "<<py[i]<<"  "<<yprime[i]<<'\n';
  istart++; to=t+tinc;
  while(to<term){
     ddassl_(isores,NN,t,py,yprime,to,info,rtol,atol,idid, rwork, lrw, iwork, liw, rpar, ipar, jac);
     if(idid==-1) { atol=1.0e-7; distr(py, yprime); info[0]=1; idid=0;
       ddassl_(isores,NN,t,py,yprime,to,info,rtol,atol,idid, rwork, lrw, iwork, liw, rpar, ipar, jac); }
     if(idid<-1) { throw("dassl problem");}
     sdev[istart] = nadh; fiout(to,fkin,istart-1); 
     t=to;
     istart++; to += tinc;
     }
  write("i0");
  distr(py, yprime);
  double fm=0;
  for(int i=0;i<NN;i++) if(fabs(yprime[i])>fabs(fm)) fm=yprime[i];
  cout<<"fm="<<fm<<'\n';
     return istart;}

void Ldistr::chast(double *py, double tint) {
  ostringstream fkin,fkont; ifin1=0;//  conc[nc2ros]=0.; conc[nc3ros]=0.;
// py[nmet+nglu]=nv.nv[suc_o];
//  ifin1= ddisolve(0,tint,py,fkin);
//  setconc(nca,20.);
//  ifin1= ddisolve(ifin1,0.5,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,1.0,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,1.5,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,2.,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,2.5,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,3.0,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,3.5,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,4.,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,4.5,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,5.,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,5.5,py,fkin); setconc(nca,20.);
//  double s=nv.setval(vlk,15.); // leak
//  ifin1= ddisolve(ifin1,2.0,py,fkin);
//  double a7=nv.setval(ptp, 0.1);	//PTP open 
//	ifin1= ddisolve(ifin1,3.0,py,fkin);
//  	nv.setval(ptp, a7);		//PTP close
//  double a6=nv.setval(qHbnd, 0.0000);//myxothiazol
//  ifin1= ddisolve(ifin1,3.0,py,fkin);
//  double a5=nv.setval(vgluout, 0.); // glutamate
//  ifin1= ddisolve(ifin1,100.0,py,fkin);
//  nv.setval(glu_o, 0.04);
//  double a1=nv.setval(qnbnd, 0.0);//antimycine
//  ifin1=ddisolve(ifin1,4.0,py,fkin);
//  double a6=nv.setval(qHbnd, 0.0000);//myxothiazol
	ifin1= ddisolve(ifin1,tint,py,fkin);
//	tisolve(tint,py,99,fkin);
//	fiout(tint,fkont,1); 
cout <<" qh= "<< conc[nqh]<<" psi= "<< conc[npsi]<<" nad= "<< conc[nnad]; 
	cout <<" pyr= "<< conc[npyr]<<"\tlac= "<< conc[nlac]<<"\n";// 
	cout <<" suc= "<< conc[nsuc]<<"\tfum= "<< conc[nfum]<<"\n";// 
	cout <<" AKG= "<< conc[nakgm]<<"\tOAA= "<< conc[noaa]<<"\n";// 
	cout <<" glu= " << conc[nglu]<<"\tgluout= " << conc[ngluo]<<"\n"; 
	cout <<" Asp= " << conc[naspm]<<"\tQ= " << conc[nqh]<<"\n"; 
//  nv.setval(vlk,s); // leak
//  nv.setval(vgluout, a5); // glutamate
//  nv.setval(qnbnd, a1);//antimycine
//  nv.setval(qHbnd, a6);//myxothiazol
//	fkin.str(std::string());
//	fkont.str(std::string());
//	ifin1= ddisolve(0,tint,py,fkin);
	fiout(tint,fkont,1); 
	write("ii");
	kont = fkont.str();
  kkin=fkin.str();
  }

inline void zerowsmal(Vec_DP& w, double cond, int n) {
	double wmax=0.0; int i;
	for (i=0;i<n;i++) if (w[i] > wmax) wmax=w[i];
	double wmin=wmax*cond;
// zero the "small" singular values
	for(i=0;i<n;i++) if(w[i] < wmin) w[i]=0.0;
}
inline double fndmax(double *ar,int n){
	double cmax(0.);
	for (int i=0;i<n;i++) {
	    if(fabs(ar[i])>cmax) cmax=fabs(ar[i]);
	    if(cmax>1e4) {cout<<"!error: max rhs="<<cmax <<'\n';
	          throw "error";}
	  }
	return cmax;
}
inline void matcp(Mat_DP& acp,Mat_DP& asr,int n){
	int i,j;
	for(i=0;i<=NN;i++)
	   for(j=0;j<=NN;j++) acp[i][j]=asr[i][j];
}
double matrext(int ipar,double dp,double ds,double *y,double *dy,double *y0,double *rhs0,double **df,double p0){
	double sum(0.),rhsmax;
	int i;
	for(i=0;i<=NN;i++) {y[i]+=dy[i]*ds; sum+=(y[i]-y0[i])*dy[i];} 
	horse.distr(y, rhs0);
	cout << "\n max rhs="<< (rhsmax=fndmax(rhs0,NN)) << '\n';
// last column in the matrix and rhs
	for(i=0;i<=NN;i++) {
	  df[i][NN]=rhs0[i]/dp; rhs0[i] *= -1.;}
	rhs0[NN]=sum-(horse.nv.nv[ipar]-p0)*dy[NN]-ds;
	sum=0.;
// last row in the matrix
	for(i=0;i<=NN;i++) sum += dy[i]*dy[i]; sum=sqrt(sum);
	for(i=0;i<=NN;i++)  df[NN][i]=dy[i]/sum; sum=0.;
  return rhsmax;}
double matdir(int ipar,double dp,double *y,Vec_DP& dy,Vec_DP& rhs,double **df){
    double x, sum(0.),rhsm;
    int i;
	Jacobian(x, y, df);
	horse.nv.nv[ipar]+=dp;
	horse.distr(y, &rhs[0]);
	cout << "\n max rhs0="<<(rhsm=fndmax(&rhs[0],NN)) << '\n';
// last column in the matrix and rhs
	for(i=0;i<NN;i++) {
	  df[i][NN]=rhs[i]/dp; rhs[i]=0;}
	rhs[NN]=1.;
// last row in the matrix
	for(i=0;i<NN;i++) {dy[i]/=dp; sum += dy[i]*dy[i];}
	sum=sqrt(sum+1.);
	for(i=0;i<NN;i++)  df[NN][i]=dy[i]/sum;
	df[NN][NN]=1./sum; 
   return rhsm;
}

void svdsolve(Mat_DP& u,Vec_DP& w,Mat_DP& v,Vec_DP& rhs,Vec_DP& dy){
	int i;
// SVD of u
	NR::svdcmp(u,w,v);
	zerowsmal(w, 1e-9, NN+1);
//	cout<<"\nw: ";
//	for (i=0;i<=NN;i++) cout<<w[i]<<" ";cout<<'\n';
// New direction vector	  
	NR::svbksb(u,w,v,rhs,dy); 
	cout<<"dy: ";
	for (i=0;i<=NN;i++) cout<<dy[i]<<" ";
	cout<<'\n';
}

void arch(double *y,int ipar,double pint,double step){
   Vec_DP w(NN+1), rhs(NN+1), dy(NN+1);
   Mat_DP a(NN+1,NN+1), u(NN+1,NN+1),v(NN+1,NN+1);
   double x, dp, y0[NN+1], *df[NN+1], ddp(1.), rhm, rhm0, sum(0);
   int i,j;
	dp=step;
// initiation: parameter continuation
   init_ydot(y,ipar,dp,dy,&rhs[0],NN);
	cout<<"*** pass init_ydot***";
// matrix for direction vector
	for(i=0;i<=NN;i++) df[i]=&u[i][0];
	double p0=horse.nv.nv[ipar],ds(0);
	Jacobian(x, y, df);
	horse.nv.nv[ipar]+=dp;
	horse.distr(y, &rhs[0]);
	cout << "\n max rhs0="<<(rhm0=fndmax(&rhs[0],NN)) << '\n';
// last column in the matrix and rhs
	for(i=0;i<NN;i++) {
	  df[i][NN]=rhs[i]/dp; rhs[i]=0;}
	rhs[NN]=1.;
// last row in the matrix
	for(i=0;i<NN;i++) {dy[i]/=dp; sum += dy[i]*dy[i];}
	sum=sqrt(sum+1.);
	for(i=0;i<NN;i++)  df[NN][i]=dy[i]/sum;
	df[NN][NN]=1./sum; 

	matcp(a,u,NN);// copy u
	svdsolve(u, w, v, rhs, dy);
	matcp(u,a,NN);// copy a
// Neuton method
// initial values for archlength continuation
	for(i=0;i<NN;i++) y0[i]=y[i];
	y0[NN]=y[NN]=p0; horse.nv.nv[ipar]=p0;
	ds=0.0001;
	for(i=0;i<9;i++){
	  dp=dy[NN]*ds;
	  horse.nv.nv[ipar]+=dp;
	  rhm=matrext(ipar,dp,ds, y, &dy[0], y0, &rhs[0], df, p0);
	  if(rhm>rhm0) {horse.nv.nv[ipar]-=dp;
	      for(i=0;i<=NN;i++) y[i]-=dy[i]*ds; break;}
	  else {svdsolve(u, w, v, rhs, dy); rhm0=rhm;}
//	  matcp(u,a,NN);// copy a
	}
	cout<<"end 1st correction"<<'\n';
	for(i=0;i<NN;i++) y0[i]=y[i]; p0=horse.nv.nv[ipar];
	y0[NN]=y[NN]=p0;
	for(i=0;i<NN;i++) rhs[i]=0; rhs[NN]=1;
	matcp(u,a,NN);// copy a
	svdsolve(u, w, v, rhs, dy);
	matcp(u,a,NN);// copy a
	
	ds=0.0003;
	for(i=0;i<9;i++){
	  dp=dy[NN]*ds;
	  horse.nv.nv[ipar]+=dp;
	  rhm=matrext(ipar,dp,ds, y, &dy[0], y0, &rhs[0], df, p0);
	  if(rhm>(rhm0+1e-3)) {horse.nv.nv[ipar]-=dp;
	      for(i=0;i<=NN;i++) y[i]-=dy[i]*ds;
	      horse.distr(y, &rhs[0]);
	      cout<<"increase! rhm="<<(rhm0=fndmax(&rhs[0],NN))<<endl;
	      break;}
	  else {svdsolve(u, w, v, rhs, dy); rhm0=rhm;}
//	  matcp(u,a,NN);// copy a
	}
}

void init_ydot(double *y,int ipar,double dp,Vec_DP& dy,double *c1,const int n){
	Vec_DP c(n), w(n);
	Mat_DP u(n,n),v(n,n);
	int i,j;
	double x, sum(0), y0[n];
	correct(y,c, w, u, v,n);
	for (i=0;i<n;i++) y0[i]= y[i]; //{c0[i]=(c1[i]-c0[i])/dp;}
	double p0=horse.nv.nv[ipar];
// last column c1/dp, rhs= -c1 :
	horse.nv.setval(ipar,p0+dp);
// last row ydot;
	correct(y,c, w, u, v,n);
	for(i=0;i<n;i++) dy[i]=(y[i]-y0[i])/dp;
	cout<<"ydot: ";
	for(i=0;i<n;i++) cout<<dy[i]<<" "; cout<<'\n';
}
void correct(double *y,Vec_DP& rhs,Vec_DP& w,Mat_DP& u,Mat_DP& v,const int n){
	int i,k;
	Mat_DP a(n,n);
	jsvd(y, w, u, v, a,1e-11,n);
//	return;
	cout<<"---correct SVD: w=\n";
	for (i=0;i<n;i++) cout<<w[i]<<" "; cout<<'\n';
	for(k=0;k<15;k++){ horse.distr(y, &rhs[0]);
	  for(i=0;i<n;i++) rhs[i] *= (-1.);
	  cout <<k;
	  if(neuton(y,rhs, w, u, v, n)<1e-7) break;
	}
}
void jsvd(double *y,Vec_DP& w,Mat_DP& u,Mat_DP& v,Mat_DP& a, double cond,const int n){// get jacobian u, copy it into a, make SVD of it as u, w, v, zero small w values:

	double x, *df[n];
	int i,j;
        for(i=0;i<n;i++) df[i]=&u[i][0];
	Jacobian(x,y,df);
// copy jacobian into a
	for (i=0;i<n;i++) //{ cout<<()<<" "; cout<<endl;}
	  for (j=0;j<n;j++) a[i][j]= u[i][j];
	NR::svdcmp(u,w,v);
	for (i=0;i<n;i++) cout<<w[i]<<" "; cout<<'\n';
	zerowsmal(w, cond, n);
//// check: a=u×w×v*
//	for(j=0;j<n;j++) for(i=0;i<n;i++) u[j][i] *= w[i];
//	for(i=0;i<n;i++) for(j=0;j<n;j++) {a[i][j]=0.;
//	   for(int k=0;k<n;k++) a[i][j]+=u[i][k] * v[j][k];}
//	cout<<"\n\n u×w×v*:\n";
//	for (i=0;i<n;i++){
//	  for (j=0;j<n;j++) cout<<(a[i][j])<<" "; cout<<endl;} 
}
double neuton(double *y,Vec_DP& c,Vec_DP& w,Mat_DP& u,Mat_DP& v,const int n){
//	cout << fixed << setprecision(6);
	Vec_DP dx(n),c1(n),dx1(n); int i,j;
	NR::svbksb(u,w,v,c,dx);
	double wmax=0.0; 
//	for ( i=0;i<n;i++){ c1[i] = 0.;
//	   for ( j=0;j<n;j++) c1[i] += a[i][j]*dx[j];
//	   c1[i] = c[i]-c1[i];}
//	cout << " max of c1: "<<fndmax(&c1[0],n) << endl;
////	jsvd(y, w, u, v, a,1e-12,n);
//	NR::svbksb(u,w,v,c1,dx1); dx[i] -= dx1[i];
	for ( i=0;i<n;i++) { y[i]+=dx[i];}
	cout << " max of dx: "<<fndmax(&dx[0],n) << endl;
	horse.distr(y, &c[0]); wmax=fndmax(&c[0],n);
	cout << " max rhs="<<wmax << endl;
return wmax; }

void isT::Function(double x, double *y, double *f){
	horse.distr(y, f);
}
void Function(double x, double *y, double *f){
	horse.distr(y, f);
}
//void Jacobian(double x, double *y, double **J){}
void Jacobian(double x, double *y, double **dfdy){
	int i,j;
//	int NN=y.size();
        double dydx0[NN];
        double dydx1[NN];
	double dy,aa;
		horse.distr(y, dydx0);
	for ( i=0;i<NN;i++) { aa=y[i];
	  dy=aa*0.001; if(dy>1.) dy=1.;
	  y[i] += dy;
	  horse.distr(y, dydx1);
	  y[i] = aa;
	  for ( j=0;j<NN;j++) 
	     if(dy>1e-12) dfdy[j][i] = (dydx1[j] - dydx0[j]) / dy;
	     else dfdy[j][i] = 0.0;
	}
}

void Mass(double **M){} // Mass
double tisolve(const double tmax,double *yy,int numint,ostringstream& fkin){
	// dimension of problem
	// initial value for x
	double xbeg(0.0), xend = tmax*tfac, dx(tmax*tfac/numint);
	// rtoler and atoler are scalars
	int itoler(0);
	// relative tolerance
	double *rtoler = new double(3.0e-7);
	// absolute tolerance
	double *atoler = new double(1.0e-12);
	const int iout(1);
	// initial step size
	double hinit(0.0005);
	// analytical Jacobian function provided
	 int ijac(0);
	// number of non-zero rows below main diagonal of Jacobian
	int mljac(NN);
	// number of non-zero rows above main diagonal of Jacobian
	int mujac(NN);
	// Mass matrix routine is identity
	const int imas(0);
	int mlmas(NN);
	int mumas(0);
	
	// Use default values (see header files) for these parameters:
	double hmax(0.0);
	int nmax(0);
	double uround(0.0), safe(0.), facl(0.0), facr(0.0);
	int nit(0);
	bool startn(false);
	int nind1(0), nind2(0), nind3(0), npred(0), m1(0), m2(0);
	bool hess(false);
	double fnewt(0.0), quot1(0.0), quot2(0.0), thet(0.1);
	
	double beta = 0.007;
	int nstiff = -1;
	int nrdens = NN;
	unsigned *icont = NULL;
	int i;
	for (i=0;i<numint;i++) {
	xend=xbeg+dx;
	isT stiffT(NN, yy, xbeg, xend, dx, itoler, rtoler, atoler,
		iout, hinit, hmax, nmax, uround, safe, facl, facr, ijac, mljac,
		mujac, imas, mlmas, mumas, nit, startn, nind1, nind2, nind3, npred,
		m1, m2, hess, fnewt, quot1, quot2, thet);
//	cout << "\n\n*******Problem integrated with RADAU5*******\n\n";
	stiffT.Integrate();
     horse.fiout(xend,fkin);
	xbeg=xend;
	}

	
	delete rtoler;
	delete atoler;
	return 1.;
}
void derivsl(const DP x, Vec_IO_DP &y, Vec_O_DP &dydx){
	DP *py=&y[0]; DP *pdydt=&dydx[0];
        nrhs++;
	horse.distr(py, pdydt);
}

double integrbs(const double tmax,Vec_DP &ystart,const int KMAX){
	int i,j; char fn[11]="init";
        DP eps=1.0e-6,h1=0.0001,hmin=1.0e-09,x1=0.0;
        nrhs=0;
        dxsav=tmax/(KMAX-1);
        kmax=KMAX;
        xp_p=new Vec_DP(KMAX);
        yp_p=new Mat_DP(NN,KMAX);
        Mat_DP &yp=*yp_p;
	Vec_DP &xp=*xp_p;
        int nbad,nok;
//	double *py = &ystart[0];
//	horse.setyinit(py);
       cout<<"pass*****"<<endl;
     try {
        NR::odeint(ystart,x1,tmax,eps,h1,hmin,nok,nbad,derivsl,NR::rkqs);
    }catch( char const* str ){ cout << "exception: "<< str <<endl; 
//    			horse.kinet(xp,yp,kount);
                         throw str;
                         }
//    else if(method==1) NR::odeint(ystart,x1,tmax,eps,h1,hmin,nok,nbad,derivsl,NR::bsstep);if(method==0)
//    else if(method==2) NR::odeint(ystart,x1,tmax,eps,h1,hmin,nok,nbad,derivsl,NR::stifbs);
//    else NR::odeint(ystart,x1,tmax,eps,h1,hmin,nok,nbad,derivsl,NR::stiff);
//    horse.label(ystart,0.1);
//    horse.kinet(xp,yp,kount);
//    for (int i=0;i<NN;i++) yy[i]=yp[i][5];
//    DP *pyy = &yy[0];
//    horse.gety(pyy);
//    horse.label();
        delete yp_p;
        delete xp_p;
return tmax;}
/*
void sens(Vec_DP &ystart,float factor){
   double *py = &ystart[0];
      	horse.setisot(py); 
    ofstream fi("total"); int ii=0;
	horse.tout(0,ii,fi);
  for (int j=1;j<horse.nv.par[0];j++){
    cout<<"parameter "<<horse.nv.par[j]<<" "<<horse.nv.pname[horse.nv.par[j]]<<endl;
   for (int i=0;i<1;i++) {
   horse.nv.read();
    horse.nv.changeval(horse.nv.par[j],factor);
     try {
       solve(py);
    }catch( char const* str ){ cout << "exception: "<< str <<endl; 
                         break; }
	horse.tout(horse.nv.par[j],ii,fi);
	} }
}*/

