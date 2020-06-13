#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
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
  double t=tex[istart]*tfac,to,rtol=0.001,atol=1.0e-9,h0=1e-5, hmax=2., rpar[2], yprime[NN], rwork[200000];
   rwork[1]=hmax; rwork[2]=h0;
     for(int i=0;i<15;i++) info[i]=0;   info[6]=1; info[10]=1;
   distr(py, yprime);
//   for(int i=0;i<NN;i++) cout<<i<<") "<<py[i]<<"  "<<yprime[i]<<'\n';
   istart++; to=tex[istart];
  while(to<tfin){
   to *= tfac;
 ddassl_(isores,NN,t,py,yprime,to,info,rtol,atol,idid,rwork,lrw,iwork, liw, rpar, ipar, jac);
 if(idid==-1) { atol=1.0e-7; distr(py, yprime); info[0]=1; idid=0;
  ddassl_(isores,NN,t,py,yprime,to,info,rtol,atol,idid,rwork,lrw,iwork, liw, rpar, ipar, jac); }
   if(idid<-1) { throw("dassl problem");}
    sdev[istart] = nadh; fiout(to,fkin,istart-1); 
     t=to;// cout<<to/tfac<<"; fc1="<<nv.flx[fc1]*tfac/60./0.4<<endl;
   istart++; to=tex[istart]; // cout<<cIIq.total(14)<<" "<<nv.frw[vqbind]/nv.nv[rct]<<endl;
 }
       write("i0");
     return istart;}

void Ldistr::chast(double *py, double tint) {
  ostringstream fkin,fkont; conc[nc2ros]=0.; conc[nc3ros]=0.;
  ifin1= ddisolve(0,tint,py,fkin);
//  setconc(nca,20.);
//  ifin1= ddisolve(ifin1,1.5,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,2.,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,2.5,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,3.0,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,3.5,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,4.,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,4.5,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,5.,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,5.5,py,fkin); setconc(nca,20.);
//  ifin1= ddisolve(ifin1,tint,py,fkin);
//  double s=nv.setval(vlk,15.);
//      double a5=nv.setval(vfadf, 0.15); // glutamate
//  ifin2= ddisolve(ifin1,2.0,py,fkin);
//      double a1=nv.setval(qnbnd, 0.0);//antimycine
//  double ifin3=ddisolve(ifin2,3.5,py,fkin);
//      double a6=nv.setval(qHbnd, 0.0000);//myxothiazol
//      conc[nc2ros]=0.; conc[nc3ros]=0.;
//      ifin1= ddisolve(0,tint,py,fkin);
      fiout(tint,fkont,1); 
////  nv.setval(vlk,s);
//  nv.setval(vfadf, a5);
//  nv.setval(qnbnd, a1);
//  nv.setval(qHbnd, a6);//myxothiazol
  kont = fkont.str();
  kkin=fkin.str();}

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
	   if(aa>0.0001) dy=aa*0.01; else dy=0.000001;
		y[i] += dy;
		horse.distr(y, dydx1);
		y[i] = aa;
		for ( j=0;j<NN;j++) {
		dfdy[j][i] = (dydx1[j] - dydx0[j]) / dy;
		}
	}
}

void Mass(double **M){} // Mass
double tisolve(const double tmax,double *yy,int numint){
	// dimension of problem
	// initial value for x
	double xbeg(0.0), xend = tmax, dx(tmax/numint-1);
	// rtoler and atoler are scalars
	int itoler(0);
	// relative tolerance
	double *rtoler = new double(3.0e-2);
	// absolute tolerance
	double *atoler = new double(1.0e-5);
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

	isT stiffT(NN, yy, xbeg, xend, dx, itoler, rtoler, atoler,
		iout, hinit, hmax, nmax, uround, safe, facl, facr, ijac, mljac,
		mujac, imas, mlmas, mumas, nit, startn, nind1, nind2, nind3, npred,
		m1, m2, hess, fnewt, quot1, quot2, thet);

		//	cout << "\n\n*******Problem integrated with RADAU5*******\n\n";
 horse.read();
	stiffT.Integrate();
//	ofstream fi("kinetics");
//     horse.fiout(xend,fi);
	
	delete rtoler;
	delete atoler;
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

