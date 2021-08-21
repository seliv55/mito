//---------------------------------------------------------------------------
#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "modlab.h"
//---------------------------------------------------------------------------
using namespace std;
using namespace Label;
double k_i= 144.53, Na_i= 0.01875, psio=86.266;
double Ldistr::o2deliv(double Pm){
  double vo2=350.,x0=20.,km=20.,v=0.83,vb=2.,kt=0.15,kl=0.189, kv=5.482, pat=160., Cmax=203.;
 pal = pat - vo2/kv,
 pa = ((-1 + kl)*kv*pat + kl*(-1 + kt)*kv*Pm + vo2 - kl*vo2)/((-1 + kl*kt)*kv), 
 pv = ((-kv)*Pm + kt*(kv*((-1 + kl)*pat + Pm) + vo2 - kl*vo2))/((-1 + kl*kt)*kv),
 vo2b = Cmax*km*(-1 + kl*kt)*kv*v*vb*(1/(km*(-1 + kl*kt)*kv + kv*(-Pm + x0) + kt*(vo2 - kl*vo2 + kv*((-1 + kl)*pat + Pm - kl*x0))) + 1/(km*(kv - kl*kt*kv) + (-1 + kl)*vo2 + kv*(pat - kl*pat - x0 + kl*(Pm - kt*Pm + kt*x0))));
return vo2b;
}
void Ldistr::seteq( double *py,double *pdydt) {
//	conc = &py[nmet];
//	dconc = &pdydt[nmet];
	setisot(py);
	setdisot(pdydt);
//  pBC1q.total(); qhpBC1.total(); BC1qn.total(); BC1.total(); cIq.total();
//  coreI.total(); cIIq.total(); coreII.total(); 
//  bc15=1.-pBC1q.cont-qhpBC1.cont-BC1qn.cont-BC1.cont;
//  coreI.e6 = 1.-coreI.cont-cIq.cont;
//  qq=qt-2*(pBC1q.cont)*c3t-qhpBC1.cont*c3t-BC1qn.cont*c3t-conc[nqh]-cIq.cont*c1t-cIIq.cont*c2t;
//  cr11 = 1.-coreII.cont-cIIq.cont;
	qq=3.75-conc[nqh];
	for (int i=0;i<NN;i++) pdydt[i]=0.;
	hfi = hi*exp(-0.255*frt*conc[npsi]);
	hfo = ho*exp(0.255*frt*conc[npsi]);
	nadh= 1. - conc[nnad];
	adp= tan - conc[nAtp];
	 }
	 
void Ldistr::c3calc( double *py,double *pdydt) {
  double kr= nv.frw[rqp_FS]*hfo*hfo;
  double ox=0.955, bp=0.;
//Qo->FeS
double sum= pBC1q.shiftFeS(nv.frw[qp_FS],kr,nv.frw[qp_bl],nv.frw[rqp_bl], nv.frw[vros], ox, dconc[nc3ros],c3t);
sum += qhpBC1.shiftFeS(nv.frw[qp_FS],kr,nv.frw[qp_bl],nv.frw[rqp_bl], nv.frw[vros], ox, dconc[nc3ros],c3t);
//bp = pBC1q.bypass(nv.frw[bypas]);
//bp += qhpBC1.bypass(nv.frw[bypas]); sum -= bp;
//nv.flx[fros]=cros;
nv.flx[fbp]=bp;
		dconc[npsi] += sum*fc*c3t;
//1: e-transport from FeS-protein (positions 3 or 1) to c1:
   sum=pBC1q.shift1(3,nv.frw[FS_c1],nv.frw[rFS_c1]);
   sum += qhpBC1.shift1(3,nv.frw[FS_c1],nv.frw[rFS_c1]);
   sum += BC1.shift1(1,nv.frw[FS_c1],nv.frw[rFS_c1]);
   sum += BC1qn.shift1(1,nv.frw[FS_c1],nv.frw[rFS_c1]);
		dconc[npsi] += sum*fc*c3t;

//3: e-transport between the hemes of cyt.b //bL->bH
 double a(exp(-0.55*frt*conc[npsi])), b(exp(0.55*frt*conc[npsi]));
 double kf= nv.frw[bl_bh] * a;
	kr= nv.frw[rbl_bh] * b;
	pBC1q.shift1(5,kf,kr);
	qhpBC1.shift1(5,kf,kr);
	BC1.shift1(3,kf,kr);
	BC1qn.shift1(3,kf,kr);
	      
//4, 5: e-transport further to q on n-side
  kf = nv.frw[bh_qn2] * hfi * hfi;
  sum = pBC1q.bhqn(6,nv.frw[bh_qn1],nv.frw[rbh_qn1],kf,nv.frw[rbh_qn2]);
  sum += BC1qn.bhqn(4,nv.frw[bh_qn1],nv.frw[rbh_qn1],kf,nv.frw[rbh_qn2]);
//		dconc[nhi] -= 2.*sum/buf/vol; 
  dconc[npsi] += sum*fc*c3t;
		
//6: binding qh2 on the p-side
  dconc[nqh] += BC1.qphbind(qhpBC1,conc[nqh],nv.frw[qHbnd],nv.frw[rqHbnd],bc15,1) *c3t;
  dconc[nqh] += BC1qn.qphbind(pBC1q,conc[nqh],nv.frw[qHbnd],nv.frw[rqHbnd])*c3t; 

//9: dissociation of qh2 on the n-side
  dconc[nqh] += BC1qn.dissqn(BC1,conc[nqh],nv.frw[qhnds],nv.frw[rqhnds],bc15,1) *c3t;
  dconc[nqh] += pBC1q.dissqn(qhpBC1,conc[nqh],nv.frw[qhnds],nv.frw[rqhnds])*c3t;

//7: binding of q on the n-side	
    BC1.bindqn(BC1qn,qq,nv.frw[qnbnd],nv.frw[rqnbnd],bc15,1);    
    qhpBC1.bindqn(pBC1q,qq,nv.frw[qnbnd],nv.frw[rqnbnd]);

//8: dissociation of q on the p-side 	
 qhpBC1.qpdiss(BC1,qq,nv.frw[qpdis],nv.frw[rqpdis]);
 pBC1q.qpdiss(BC1qn,qq,nv.frw[qpdis],nv.frw[rqpdis]);

//12: e-transport out from c1 (assumed to c)
        kf = nv.frw[vc1c]*a;
	sum = pBC1q.c1c(4,kf, ox); 
	sum += BC1.c1c(2,kf, ox,1);
	sum += BC1qn.c1c(2,kf, ox);
	sum += qhpBC1.c1c(4,kf, ox);
	nv.flx[fc1c] = sum;
    dconc[npsi] += 2.*sum*fc*c3t;
    }
    
void Ldistr::tca(double *py,double *pdydt) {
  double v, v1, vfac(1.);
// TCA cycle from oaa to akg
	v =nv.flx[fcs] = nv.frw[vcs]*conc[noaa]*conc[npyr]*conc[nnad];
	v1=v*tnad;
	dconc[noaa] -= v1/vfac;    
	dconc[npyr] -= v1/vfac;  
	dconc[nakgm]+= v1/vfac;  
	dconc[nnad] -= 2.*v;

// Malic enzyme fum->pyr
	v = nv.frw[vmalic]*conc[nnad]*conc[nfum];
	dconc[nfum] -= v;    
	dconc[npyr] += v;  
	dconc[nnad] -= v;

// akg to succinate
	v = nv.frw[akgsuc]*conc[nnad]*conc[nakgm];
	v1=v*tnad;
	dconc[nakgm]-= v1/vfac;   
	dconc[nnad] -= v;  
	dconc[nsuc] += v1/vfac;  

// glutamate to akg (Glutamate dehidrogenase):
	v = nv.frw[vglu_akg]*conc[nglu]*conc[nnad];
	dconc[nakgm] += v/vfac;/**/
	dconc[nglu] -= v/vfac;  
	dconc[nnad] -= v;  
//leak
	v = nv.flx[flk] = 0.005*conc[npsi]*leak(nv.frw[vlk]);
	dconc[npsi] -= 2.*v*fc;
// complex I
	v = nv.frw[vcI]*nadh*tnad*qq* exp(-0.2*frt*conc[npsi]);
	dconc[nnad] += v;
	dconc[nqh]  += v;
	dconc[npsi] += 8.*v*fc;
// succinate dehydrogenase, complex II
	v = nv.flx[fsdh] = nv.frw[vsdh]*qq*conc[nsuc]*(1/(1+conc[noaa]/0.0002));
	dconc[nqh] += v; 
	dconc[nsuc] -= v; 
	dconc[nfum] += v;
// complex III
	v = nv.frw[vcIII]*qq*conc[nqh]* exp(-0.2*frt*conc[npsi]);
	dconc[nqh] -= v; 
	dconc[npsi] += 4.*v*fc;
	nv.flx[fc1c] = v;
	
//   if(conc[nca]> 0.001) {
//     double a=conc[npsi]*frt;
//     double e=exp(a);
//     v = 0.0005*a*conc[nca]*e/(1-e);
//     dconc[nca] += v;
//     dconc[npsi] += 4.*v*fc;
//   }
    }

void Ldistr::shutl(){ 
// MDH
  double x = nv.frw[vmdh]*(conc[nnad]*conc[nfum]-nadh*conc[noaa]);//
  double x1=x*tnad;
      dconc[nnad] -= x; 
      dconc[nfum] -= x1; 
      dconc[noaa] += x1; 
// aspartat aminotransferase: oaa + glum <-> aspm + akg
	x = nv.frw[vasp_atf]*conc[noaa]*conc[nglu] - nv.frw[vasp_atr]*conc[naspm]*conc[nakgm];
	dconc[noaa] -= x;
	dconc[nglu] -= x;
	dconc[naspm] += x;
	dconc[nakgm] += x;
// aspartat efflux asp ->  
	x= nv.frw[vaspout]*conc[naspm];
	dconc[naspm] -= x;
}

void Ldistr::glycolysis(){
   // glucose to pyruvate
   double x= nv.frw[vGl]*(1-conc[npyr])*conc[nnad];
   double x1=x*tnad;
     dconc[npyr] += x1;
     dconc[nnad] -= x;
   // pyruvate to lactate
   x= nv.frw[vldh]*(conc[npyr]*nadh - conc[nlac]*conc[nnad]);
   x1=x*tnad;
      dconc[npyr]  -= x1;
      dconc[nlac] += x1;
      dconc[nnad]  += x;
   // lactate efflux/uptake
   x= nv.frw[vlacdif]*(nv.nv[laco] - conc[nlac]);
      dconc[nlac] += x;
}

void Ldistr::ions(int nci,double P){ //ion equilibrium
   double mu=psio*frt;
   double e=exp(mu);
   double x=mu*P*(k_i-nv.nv[k_o]*e)/(1-e);//0;//
//   dconc[nci] += x;
}

void Ldistr::glufl(){ // equilibrium flux of external glutamate transport;
   const int len=500; double yy[500], x[len], y[len], mm(0.);
   setdisot(yy);
   stringstream gluT;
   conc[ngluo]=0.; psio=80;
   for (int i=0;i<len;i++){
      for (int j=0;j<1000;j++) jglu0();
      x[i]=conc[ngluo]; y[i]=jglu0();
      cout<<x[i] <<"  "<<y[i]<<'\n';
      conc[ngluo] +=0.00001;
   }
   for (int i=0;i<len;i++) if(y[i]>mm) mm=y[i];
   for (int i=0;i<len;i++) y[i] = y[i]/mm;
   for (int i=0;i<len;i++) gluT<<x[i]<<"  "<<y[i]<<'\n';
   ofstream fglu("glutr");
   fglu<<gluT.str();
}

inline double N3ToHG(double To,double Glu,double Na, double K){
   double a=125*Glu*Na*Na*Na;
   return (a*To)/(a+(12500*Glu+572)*Na*Na+5500*Na+110*K+1375000);
}

inline double TiK(double Ti,double Glu,double Na, double K){
   double a=50*Glu*Na*Na*Na;
   return (103180*K*Ti)/(a+(2340*Glu+10318)*Na*Na+25795*Na+103180*K+103180);
}

double Ldistr::jglu0(){ // transport of external glutamate
   double K1Na=50, K2Na=8.4, Kglu=0.0025, kt0=nv.frw[vglu_tr]*(1/(1+conc[nglu]/0.1)), kr0=nv.frw[vglu_tr]*(1/(1+conc[nglu]/0.1)), Tg;
//   Kd = Kglu*(K1Na/nv.nv[Na_o] + 1.)*K2Na/(K2Na + nv.nv[Na_o]);
   Tg=N3ToHG(conc[neo],conc[ngluo],nv.nv[Na_o],nv.nv[k_o]);
   double kt=kt0*exp(2*psio*frt/2), kr=kr0*exp(-psio*frt/2);
   double x=kt*Tg;
   dconc[neo] -= x; 
   dconc[ngluo] -= x;
   dconc[nglu] += x;
//   dconc[nNa_i] += 3*x;
   double xr=kr*TiK((1.-conc[neo]),conc[nglu],Na_i,k_i);
   dconc[neo] += xr; //conc[neo] += xr/1000; conc[nei] -= xr/1000;
//   dconc[nk_i] -= xr;
   return x;
}

void Ldistr::atpsyn(double katp){
  double v=katp*adp*conc[npsi];
  dconc[nAtp] += v;
  dconc[npsi] -= 8.*v*fc;
}

void Ldistr::NaKatpase(double k){
   double v=k*conc[nAtp]*Na_i*nv.nv[k_o];
   dconc[nAtp] -= v;
//   dconc[nk_i] += 2*v;
//    dconc[nNa_i] -= 3*v;
}

void Ldistr::atpase(double k){
   double v=k*conc[nAtp];
   dconc[nAtp] -= v;
}

void Ldistr::peroxidase(double k){
   double v3=k*nadh*conc[nc3ros];
   dconc[nc3ros] -= v3;
   double v2=k*nadh*conc[nc2ros];
   dconc[nc2ros] -= v2;
   double v1=k*nadh*conc[nc1ros];
   dconc[nc1ros] -= v1;
   dconc[nnad] += v1+v2+v3;
}

void Ldistr::transition(double kptp){
   double v;
//   v=kptp*conc[nakgm];
//   dconc[nakgm] -= v;
//   v=kptp*conc[noaa];
//   dconc[noaa] -= v;
//   v=kptp*conc[nfum];
//   dconc[nfum] -= v;
//   v=kptp*conc[nsuc];
//   dconc[nsuc] -= v;
   v=kptp*conc[nglu]*0.35;
   dconc[nglu] -= v;
   v = 0.005*conc[npsi]*leak(1.);
   dconc[npsi] -= 2.*v*fc;
}

void Ldistr::c1calc( double *py,double *pdydt) {
//COMPLEX I (qp-qp-qn-qn-n2-fmn-fmn)
 nv.flx[fc1] = coreI.fmnred(nv.frw[vfred], nv.frw[vrfred], nadh*tnad, conc[nnad]*tnad, fmn, fmnh);
 nv.flx[fc1] += cIq.fmnred(nv.frw[vfred], nv.frw[vrfred], nadh*tnad, conc[nnad]*tnad, fmn, fmnh); 
 dconc[nnad] += nv.flx[fc1]*c1t/tnad; //NADH->FMN
 coreI.n562(nv.frw[vn56],nv.frw[vrn56],n5red, n6ar);
 cIq.n562(nv.frw[vn56],nv.frw[vrn56],n5red, n6ar);
	double kf1 = nv.frw[vn2qn1] * exp(-0.5*frt*conc[npsi]);
	double kr1 = nv.frw[vrn2qn1] * exp(0.5*frt*conc[npsi]); 
	double kf2 = nv.frw[vn2qn2];
	double kr2 = nv.frw[vrn2qn2]; 
 double sum = cIq.n2q(kf1,kr1,kf2,kr2,n2red)*c1t;// 0.;// n2 -> Q
 dconc[npsi] += 8.*sum*fc; nv.flx[fc1] = sum;
 fsq = coreI.getfs(fs) + cIq.getfs(fs) + cIq.getsq();
 dconc[nc1ros] = fsq*nv.frw[vros];
 dconc[nqh] += cIq.qhdiss1(coreI,conc[nqh],nv.frw[vndis],nv.frw[vrndis])*c1t;
 coreI.qbind1(cIq,qq,nv.frw[vpbind],nv.frw[vrpbind]);
}
void Ldistr::c2calc( double *py,double *pdydt) {
//COMPLEX II (q-q-fs-fs-fs-fad-fad)
 double k=nv.frw[vfadf]*(1/(1+conc[noaa]/0.0002));
 nv.flx[fc2] = coreII.fadred(k, nv.frw[vfadr], conc[nsuc], conc[nfum],cr11);
 nv.flx[fc2] += cIIq.fadred(k, nv.frw[vfadr], conc[nsuc], conc[nfum],cr11); 
 dconc[nfum] += nv.flx[fc2]*c2t/1; 
 dconc[nsuc] -= nv.flx[fc2]*c2t/10; 
	double kf1 = nv.frw[vbq1] ;
	double kr1 = nv.frw[vrbq1]; 
	double kf2 = nv.frw[vbq2];
	double kr2 = nv.frw[vrbq2];
 coreII.fadfs(kf1,kr1,f2s);
 cIIq.fadfs(kf1,kr1,f2s);
 double sum = cIIq.fsq(kf2,kr2,kf2,kr2,br);
 dconc[nqh] += cIIq.qhdiss1(coreII,conc[nqh],nv.frw[vqdis],nv.frw[vrqdis], cr11)*c2t;
 coreII.qbind1(cIIq,qq,nv.frw[vqbind],nv.frw[vrqbind],cr11);
 dconc[nc2ros]=coreII.fadros(nv.frw[vros])*c2t;
 dconc[nc2ros]+=cIIq.fadros(nv.frw[vros])*c2t;
// dconc[nc2ros]+=cIIq.sqros(nv.frw[c2ros]);
}

//void Ldistr::distr( double *py,double *pdydt) {
//  seteq(py,pdydt); 
//  tca(py,pdydt);
////  c3calc(py,pdydt);
////  c1calc(py,pdydt);
////  c2calc(py,pdydt);
////  ions(nk_i,vK);
//  jglu0();
//  dconc[ngluo]+=gluout(nv.frw[vgluout],nv.nv[glu_o]);
//  atpsyn(nv.frw[vatsyn]);
////  NaKatpase(nv.frw[katase]);
//  shutl();
//  glycolysis();
//  atpase(nv.frw[vatpase]);
////  peroxidase(nv.frw[vperos]);
////  if(conc[nqh]>3.3) transition(nv.frw[ptp]);
////  dconc[nsuc]=0;
////  dconc[nglu]=0;
////  dconc[nakgm]=0;
////  dconc[nfum]=0;
////  dconc[npyr]=0;
////  dconc[nnad]=0;
////  dconc[nNa_i]=0;
////  dconc[nk_i]=0;
//}
void Ldistr::distr( double *py,double *pdydt) {
double xk_i= 144.53, xNa_i= 0.01875, psio=86.266, fc=500., frt=0.039, hi=5.4e-05, ho=0.000114, krct= 1000,
xlaco = 0.3,
xglu_o = nv.nv[glu_o],
xNa_o = 140,
xk_o = 5,
klk = 0,
 katsyn = nv.frw[vatsyn],// 0.001*krct, 	// ADP -> ATP
 katpase = nv.frw[vatpase],// 1.9e-04*krct, 	// ATP -> ADP
 kGl = nv.frw[vGl],// 0.001*krct, 		// glucose input
 kldh = nv.frw[vldh],// 0.001*krct, 		// pyr<-> lac
 klacdif = nv.frw[vlacdif],// 0.0001*krct, 	// lac out <-> lac in
 kcI = nv.frw[vcI],// 0.0001*krct,		// NADH+q <-> NAD+qh
 ksdh = nv.frw[vsdh],// 0.004*krct, 		// suc+q <-> fum+qh
 kcIII = nv.frw[vcIII],// 0.004*krct, 	// qh+q <-> 2q
 kcs = nv.frw[vcs],// 0.07*krct,
 kmalic = nv.frw[vmalic],// 0.0005*krct,  	// oaa -> pyr
 kakgsuc = nv.frw[akgsuc],// 0.01*krct, 	// akg -> suc
 kmdh = nv.frw[vmdh],// 0.9*krct, 		// mal -> oaa
 kasp_atf = nv.frw[vasp_atf],// 0.45*krct, 	// oaa + glu -> asp + akg
 kasp_atr = nv.frw[vasp_atr],// 0.1*krct, 	// oaa + glu <- asp + akg
 kaspout = nv.frw[vaspout],// 0.03*krct, 		// asp efflux
//62 vK 0.001 		// ** K+ in -> K+ out
 kglu_akg = nv.frw[vglu_akg],// 0.0001*krct, 	// Glu -> akg
 kglu_tr = nv.frw[vglu_tr],// 0.3*krct, 	// Glu out -> Glu 
 kgluout = nv.frw[vgluout];// 0.075*krct; 	// spilover
//66 kATPase 0.001 	// K/Na ATPase
	setisot(py);
	setdisot(pdydt);
	qq=3.75-conc[nqh];
	for (int i=0;i<NN;i++) pdydt[i]=0.;
	hfi = hi*exp(-0.255*frt*conc[npsi]);
	hfo = ho*exp(0.255*frt*conc[npsi]);
	nadh= 1. - conc[nnad];
	adp= tan - conc[nAtp];
// TCA cycle from oaa to akg
 double v1 = kcs*conc[noaa]*conc[npyr]*conc[nnad];
// Malic enzyme fum->pyr
 double v2 = kmalic*conc[nnad]*conc[nfum];
// akg to succinate
 double v3 = kakgsuc*conc[nnad]*conc[nakgm];
// glutamate to akg (Glutamate dehidrogenase):
 double v4 = 40*kglu_akg*conc[nglu]*conc[nnad];
//leak
 double v5 = 0.005*conc[npsi]*klk*exp(frt*conc[npsi])*(ho - hi);
// complex I
 double v6 = kcI*nadh*tnad*qq* exp(-0.2*frt*conc[npsi]);
// succinate dehydrogenase, complex II
 double v7 = 10*ksdh*qq*conc[nsuc]*(1/(1+conc[noaa]/0.0002));
// complex III
 double v8 = kcIII*qq*conc[nqh]* exp(-0.2*frt*conc[npsi]);
// glutamate transport
 double a=125*conc[ngluo]*xNa_o*xNa_o*xNa_o;
 double Tg=a*conc[neo]/(a+(12500*conc[ngluo]+572)*xNa_o* xNa_o + 5500*xNa_o + 110*xk_o + 1375000);
 double kk=1;//(atan(1000*(conc[nglu]-2.5))/1.554-1)/(-2.);
 double kt0=kglu_tr*(1/(1+40*conc[nglu]/0.1));
 double kr0=kglu_tr*(1/(1+40*conc[nglu]/0.1));
 double kt=kt0*exp(2*psio*frt/2), kr=kr0*exp(-psio*frt/2);
 double v9=kk*kt*Tg;
 double nei=1.-conc[neo];
	a=50*40*conc[nglu]*xNa_i*xNa_i*xNa_i;
 double v10=kk*kr*103180*xk_i*nei/(a+(2340*40*conc[nglu]+10318) *xNa_i*xNa_i +25795*xNa_i+103180*xk_i+103180);
// ATP synthase
 double v11=katsyn*adp*conc[npsi];
// MDH
 double v12 = kmdh*(conc[nnad]*conc[nfum]-nadh*conc[noaa]);//
// aspartat aminotransferase: oaa + glum <-> aspm + akg
 double v13 = 40*kasp_atf*conc[noaa]*conc[nglu] - kasp_atr*conc[naspm]*conc[nakgm];
// aspartat efflux asp ->  
 double v14= kaspout*conc[naspm];
// glucose to pyruvate
 double v15= kGl*(1-conc[npyr])*conc[nnad];
// pyruvate to lactate
 double v16= kldh*(conc[npyr]*nadh - conc[nlac]*conc[nnad]);
// lactate efflux/uptake
 double v17= klacdif*(xlaco - conc[nlac]);
// atpase
 double v18=katpase*conc[nAtp];
 double v19=kgluout*(xglu_o-conc[ngluo]);
// PTP
 nv.frw[ptp]=0.0001*(atan(1000*(conc[nqh]-3.73))/1.57+1)/2.;
 double v20=1*nv.frw[ptp]*40*conc[nglu];
 double v21=1*nv.frw[ptp]*conc[nakgm];
 double v22=1*nv.frw[ptp]*conc[noaa];
 double v23=1*nv.frw[ptp]*conc[nfum];
 double v24=500*nv.frw[ptp]*conc[nsuc];
 double v25=1*nv.frw[vros]*conc[nqh];
	dconc[nqh] = v6 + v7 - v8 - v25;
	dconc[npsi] = 8.*v6*fc - 2.*v5*fc + 4.*v8*fc - 8.*v11*fc;
	dconc[nnad] = v6 - 2.*v1 - v2 - v3 - v4 - v12 - v15 + v16;
	dconc[npyr] = v15*tnad + v2 - v16*tnad - v1*tnad;
	dconc[nsuc] = (v3*tnad - v7 - v24)/10.;  
	dconc[nfum] = v7 - v2 - v12*tnad - v23;
	dconc[noaa] = v12*tnad - v1*tnad - v13 - v22;
	dconc[nakgm] = v1*tnad + v13 + v4 - v3*tnad - v21;
	dconc[nglu] = (v9 - v4 - v13 - v20)/40;
	dconc[ngluo] = v19 - v9;
	dconc[naspm] = v13 - v14;
	dconc[nAtp] = v11 - v18;
	dconc[neo] = v10 - v9;
	dconc[nlac] = v16*tnad + v17;
	nv.flx[fc1c] = v8;
}   

