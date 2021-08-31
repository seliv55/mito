//---------------------------------------------------------------------------
#include <iostream>
#include "nr.h"
#include "nums.hh"
#include "modlab.h"
//---------------------------------------------------------------------------
using namespace std;
using namespace Label;
double Ldistr::o2deliv(double Pm){
  double vo2=350.,x0=20.,km=20.,v=0.83,vb=2.,kt=0.15,kl=0.189, kv=5.482, pat=160., Cmax=203.;
 pal = pat - vo2/kv,
 pa = ((-1 + kl)*kv*pat + kl*(-1 + kt)*kv*Pm + vo2 - kl*vo2)/((-1 + kl*kt)*kv), 
 pv = ((-kv)*Pm + kt*(kv*((-1 + kl)*pat + Pm) + vo2 - kl*vo2))/((-1 + kl*kt)*kv),
 vo2b = Cmax*km*(-1 + kl*kt)*kv*v*vb*(1/(km*(-1 + kl*kt)*kv + kv*(-Pm + x0) + kt*(vo2 - kl*vo2 + kv*((-1 + kl)*pat + Pm - kl*x0))) + 1/(km*(kv - kl*kt*kv) + (-1 + kl)*vo2 + kv*(pat - kl*pat - x0 + kl*(Pm - kt*Pm + kt*x0))));
return vo2b;
}
void Ldistr::seteq( double *py,double *pdydt) {
	setisot(py);
	setdisot(pdydt);
    pBC1q.total(); qhpBC1.total(); BC1qn.total(); BC1.total(); cIq.total();
    coreI.total(); cIIq.total(); coreII.total(); 
  bc15=1.-pBC1q.cont-qhpBC1.cont-BC1qn.cont-BC1.cont;
  coreI.e6 = 1.-coreI.cont-cIq.cont;
  qq=qt-2*(pBC1q.cont)*c3t-qhpBC1.cont*c3t-BC1qn.cont*c3t-conc[nqh]-cIq.cont*c1t-cIIq.cont*c2t;
//  cout<<"qh2="<<conc[nqh] <<"; q="<<qq<<endl;
  cr11 = 1.-coreII.cont-cIIq.cont;
	for (int i=0;i<NN;i++) pdydt[i]=0.;
	 hfi = conc[nhi]*exp(-0.255*frt*conc[npsi]);
	 hfo = ho*exp(0.255*frt*conc[npsi]);
       nadh= 1. - conc[nnad];
       nadhc= 1. - conc[nnadc];
       adp=tan-conc[nAtp];
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
//1: e-trasport from FeS-protein (positions 3 or 1) to c1:
   sum=pBC1q.shift1(3,nv.frw[FS_c1],nv.frw[rFS_c1]);
   sum += qhpBC1.shift1(3,nv.frw[FS_c1],nv.frw[rFS_c1]);
   sum += BC1.shift1(1,nv.frw[FS_c1],nv.frw[rFS_c1]);
   sum += BC1qn.shift1(1,nv.frw[FS_c1],nv.frw[rFS_c1]);
		dconc[npsi] += sum*fc*c3t;

//3: e-trasport between the hemes of cyt.b //bL->bH
 double a(exp(-0.55*frt*conc[npsi])), b(exp(0.55*frt*conc[npsi]));
 double kf= nv.frw[bl_bh] * a;
	kr= nv.frw[rbl_bh] * b;
	pBC1q.shift1(5,kf,kr);
	qhpBC1.shift1(5,kf,kr);
	BC1.shift1(3,kf,kr);
	BC1qn.shift1(3,kf,kr);
	      
//4, 5: e-trasport further to q on n-side
  kf = nv.frw[bh_qn2] * hfi * hfi;
  sum = pBC1q.bhqn(6,nv.frw[bh_qn1],nv.frw[rbh_qn1],kf,nv.frw[rbh_qn2]);
  sum += BC1qn.bhqn(4,nv.frw[bh_qn1],nv.frw[rbh_qn1],kf,nv.frw[rbh_qn2]);
  dconc[npsi] += 2.*sum*fc*c3t;
  dconc[nhi] -= 2.*sum/buf;
		
//6: binding qh2 on the p-side
  dconc[nqh] += BC1.qphbind(qhpBC1,conc[nqh],nv.frw[qHbnd],nv.frw[rqHbnd],bc15,1)*c3t;
  dconc[nqh] += BC1qn.qphbind(pBC1q,conc[nqh],nv.frw[qHbnd],nv.frw[rqHbnd])*c3t; 

//9: dissociation of qh2 on the n-side
  dconc[nqh] += BC1qn.dissqn(BC1,conc[nqh],nv.frw[qhnds],nv.frw[rqhnds],bc15,1)*c3t;
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
//11: succinate dehydrogenase
  double v, v1, vfac(1.);
// Citrate Synthase
  v =nv.flx[fcs] = nv.frw[vcs]*conc[noaa]*conc[npyr]*conc[nnad];
  v1=v*tnad;
    dconc[noaa] -= v1/vfac;    
    dconc[npyr] -= v1/vfac;  
    dconc[ncit] += v1/vfac;  
      dconc[nnad] -= v;

//TCA cycle from cit to akg
  v =nv.flx[ftca] = nv.frw[vtca]*conc[nnad]*conc[ncit];
  v1=v*tnad;
    dconc[ncit] -= v1/vfac;   
    dconc[nnad] -= v;  
    dconc[nakgm] += v1/vfac;  

// akg to succinate
  v =nv.flx[ftca] = nv.frw[vtca]*conc[nnad]*conc[nakgm];
  v1=v*tnad;
    dconc[nakgm] -= v1/vfac;   
    dconc[nnad] -= v;  
    dconc[nsuc] += v1/vfac;  

// glutamate to akg:
  v = nv.frw[vglu_suc]*conc[nglum];
    dconc[nakgm] += v/vfac;/**/
    dconc[nglum] -= v/vfac;  
//leak
  nv.flx[flk] = 0.005*conc[npsi]*leak(nv.frw[vlk]);
  v = nv.flx[flk];
    dconc[npsi] -= 2.*v*fc;
    dconc[nhi] += v/buf;

   if(conc[nca]> 0.001) {
     double a=conc[npsi]*frt;
     double e=exp(a);
     v = 0.01*a*conc[nca]*e/(1-e);
     dconc[nca] += v;
     dconc[npsi] += 4.*v*fc;
   }
    }

void Ldistr::shutl(){ // malate-aspartate shutle
    // MDH in cytosol nadhc + oaac <-> nadc + fumc
   double x= nv.frw[vmdhc]*(nadhc*conc[noaac]);//-conc[nnadc]*conc[nfumc]
   double x1=x*tnadc;
      dconc[noaac] -= x1;
      dconc[nnadc] += x;
      dconc[nfumc] += x1;
//Fumarate oxidation and MDH reaction in mito
   x = nv.frw[vmdh]*(conc[nnad]*conc[nfum]);//-nadh*conc[noaa]
   x1=x*tnad;
      dconc[nnad] -= x; 
      dconc[nfum] -= x1; 
      dconc[noaa] += x1; 
    // transport: malate/akg antiport
    x= nv.frw[vmalakg]*(conc[nfumc]*conc[nakgm]-conc[nfum]*conc[nakgc]);
      dconc[nfumc] -= x;
      dconc[nakgm] -= x;
      dconc[nfum] += x;
      dconc[nakgc] += x;
    // transport: aspartate/glutamate antiport
    x= nv.frw[vgluasp]*(conc[nglu]*conc[naspm] - conc[nglum]*conc[naspc]);
      dconc[nglu] -= x;
      dconc[naspm] -= x;
      dconc[nglum] += x;
      dconc[naspc] += x;
      dconc[npsi] -= 2.*x*fc;
     // transport: glutamate-/OH-
     x= nv.frw[vgluOH]*(conc[nglu]-conc[nglum]);
      dconc[nglu] -= x;
      dconc[nglum] += x;
    //aspartat aminotransferase mito: oaa + glum <-> aspm + akg
    x = nv.frw[vasp_atf]*conc[noaa]*conc[nglum] - nv.frw[vasp_atr]*conc[naspm]*conc[nakgm];
      dconc[noaa] -= x;
      dconc[nglum] -= x;
      dconc[naspm] += x;
      dconc[nakgm] += x;
   //aspartat aminotransferase cytosol: oaac + glu <-> asp + akg
    x = nv.frw[vasp_atf]*conc[noaac]*conc[nglu]-nv.frw[vasp_atr]*conc[naspc]*conc[nakgc];//
      dconc[naspc]  += x;
      dconc[nakgc] += x;
      dconc[noaac] -= x;
      dconc[nglu]  -= x;
}

void Ldistr::glycolysis(){
   // glucose to pyruvate
   double x= nv.frw[vGl];
   double x1=x*tnadc;
     dconc[npyr] += x1;
     dconc[nnadc] -= x;
   // pyruvate to lactate
   x= nv.frw[vldh]*(conc[npyr]*nadhc - conc[nlac]*conc[nnadc]);
   x1=x*tnadc;
      dconc[npyr]  -= x1;
      dconc[nlac] += x1;
      dconc[nnadc]  += x;
   // lactate efflux/uptake
   x= nv.frw[vldh]*(nv.nv[laco] - conc[nlac]);
      dconc[nlac] += x;
}

void Ldistr::ions(int nci,double P){ //ion equilibrium
   double mu=conc[npsio]*frt;
   double e=exp(mu);
   double x=mu*P*(conc[nci]-nv.nv[k_o]*e)/(1-e);//0;//
   dconc[nci] += x;
   dconc[npsio] -= 2*x*fc;
}

void Ldistr::glufl(){ // equilibrium flux of external glutamate transport;
   const int len=500; double yy[500], x[len], y[len], mm(0.);
   setdisot(yy);
   stringstream gluT;
   conc[ngluo]=0.; conc[npsio]=80;
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
   double K1Na=50, K2Na=8.4, Kglu=0.0025, kt0=nv.frw[vglu_oi], kr0=nv.frw[vglu_oi], Tg, Kd;
   Kd = Kglu*(K1Na/nv.nv[Na_o] + 1.)*K2Na/(K2Na + nv.nv[Na_o]);
   Tg=N3ToHG(conc[neo],conc[ngluo],nv.nv[Na_o],nv.nv[k_o]);
   double kt=kt0*exp(2*conc[npsio]*frt/2), kr=kr0*exp(-conc[npsio]*frt/2);
   double x=kt*Tg;
   dconc[neo] -= x; 
   dconc[nei] += x; //conc[nei] += x/1000;conc[neo] -= x/1000;
   dconc[ngluo] -= x;
   dconc[nglu] += x;
   dconc[nNa_i] += 3*x;
   dconc[npsio] -= 6*x*fc;
   double xr=kr*TiK(conc[nei],conc[nglu],conc[nNa_i],conc[nk_i]);
   dconc[nei] -= xr;
   dconc[neo] += xr; //conc[neo] += xr/1000; conc[nei] -= xr/1000;
   dconc[nk_i] -= xr;
   dconc[npsio] += 2*xr*fc;
   return x;
}

void Ldistr::atpsyn(double katp){
  katp *= 1-(atan(1000*(conc[nqh]-2.93))/1.57+1)/2.;
  double v=katp*adp*conc[npsi];
  dconc[nAtp] += v;
  dconc[npsi] -= 6.*v*fc;
  dconc[nhi] += 3.9*v/buf;
}

void Ldistr::NaKatpase(double katp){
   double v=katp*conc[nAtp]*conc[nNa_i]*nv.nv[k_o];
   dconc[nAtp] -= v;
   dconc[nk_i] += 2*v;
    dconc[nNa_i] -= 3*v;
   dconc[npsio] += 2*v*fc;
}

void Ldistr::atpase(double katp){
   double v=katp*conc[nAtp];
   dconc[nAtp] -= v;
}

void Ldistr::c1calc( double *py,double *pdydt) {
//COMPLEX I (qp-qp-qn-qn-n2-fmn-fmn)
 nv.flx[fc1] = coreI.fmnred(nv.frw[vfred], nv.frw[vrfred], nadh*tnad, conc[nnad]*tnad, fmn, fmnh);
 nv.flx[fc1] += cIq.fmnred(nv.frw[vfred], nv.frw[vrfred], nadh*tnad, conc[nnad]*tnad, fmn, fmnh); 
 dconc[nnad] += nv.flx[fc1]*c1t/tnad; //NADH->FMN
 coreI.n562(nv.frw[vn56],nv.frw[vrn56],n5red, n6ar);
 cIq.n562(nv.frw[vn56],nv.frw[vrn56],n5red, n6ar);
	double kf1 = nv.frw[vn2qn1] ;
	double kr1 = nv.frw[vrn2qn1] ; 
	double kf2 = nv.frw[vn2qn2]* exp(-0.95*frt*conc[npsi])*conc[nhi]/0.0001;
	double kr2 = nv.frw[vrn2qn2]* exp(0.95*frt*conc[npsi])*ho/0.0001; 
 double sum = cIq.n2q(kf1,kr1,kf2,kr2,n2red)*c1t;// 0.;// n2 -> Q
 dconc[npsi] += 8.*sum*fc; nv.flx[fc1] = sum;
  dconc[nhi] -= 4.*sum/buf;
 dconc[nqh] += cIq.qhdiss1(coreI,conc[nqh],nv.frw[vndis],nv.frw[vrndis])*c1t; // QH2->
 coreI.qbind1(cIq,qq,nv.frw[vpbind],nv.frw[vrpbind]);// QH2<-
}
void Ldistr::c2calc( double *py,double *pdydt) {
//COMPLEX II (q-q-fs-fs-fs-fad-fad)
 double k=nv.frw[vfadf]*(1/(1+conc[noaa]/0.0002));//(1+conc[nglum]/0.013);//
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
 dconc[nc2ros]=coreII.fadros(nv.frw[c2ros])*c2t;
 dconc[nc2ros]+=cIIq.fadros(nv.frw[c2ros])*c2t;
// dconc[nc2ros]+=cIIq.sqros(nv.frw[c2ros]);
}
// PTP
void Ldistr::ptp(){
 nv.frw[vptp]=0.0025*(atan(1000*(conc[nqh]-2.93))/1.57+1)/2.;
 double v20=1*nv.frw[vptp]*40*conc[nglu]; dconc[nglu] -= v20;
 double v21=1*nv.frw[vptp]*conc[nakgm]; dconc[nakgm] -= v21;
 double v22=1*nv.frw[vptp]*conc[noaa]; dconc[noaa] -= v22;
 double v23=1*nv.frw[vptp]*conc[nfum]; dconc[nfum] -= v23;
 double v24=500*nv.frw[vptp]*conc[nsuc]; dconc[nsuc] -= v24;
 double v25=5000*nv.frw[vptp]*(ho-conc[nhi]); dconc[nhi] += v25/buf;
 dconc[npsi] -= 10000*nv.frw[vptp]*conc[npsi];
}


void Ldistr::distr( double *py,double *pdydt) {
  seteq(py,pdydt); 
  c3calc(py,pdydt);
  tca(py,pdydt);
  c1calc(py,pdydt);
  c2calc(py,pdydt);
  ions(nk_i,vK);
  jglu0();
  dconc[ngluo]+=gluout(nv.frw[vgluout],nv.nv[glu_o]);
  atpsyn(nv.frw[vatsyn]);
  NaKatpase(nv.frw[katase]);
  shutl();
  glycolysis();
  atpase(nv.frw[vatpase]);
  ptp();
//  if (conc[nhi]<1e-10) dconc[nhi]*=0.0;
//  dconc[nsuc]=0;
//  dconc[nglu]=0;
//  dconc[naspc]=0;
//  dconc[nakgc]=0;
//  dconc[noaac]=0;
//  dconc[npyr]=0;
//  dconc[nnad]=0;
}

