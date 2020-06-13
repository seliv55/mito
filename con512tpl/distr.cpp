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
    pBC1q.total(); qhpBC1.total(); BC1qn.total(); BC1.total(); cIq.total(); coreI.total(); cIIq.total(); coreII.total(); 
  bc15=c3t-pBC1q.cont-qhpBC1.cont-BC1qn.cont-BC1.cont;
  coreI.e6 = c1t-coreI.cont-cIq.cont;
  qq=qt-2*(pBC1q.cont)-qhpBC1.cont-BC1qn.cont-conc[nqh]-cIq.cont-cIIq.cont;
  cr11 = c2t-coreII.cont-cIIq.cont;
	for (int i=0;i<NN;i++) pdydt[i]=0.;
	 hfi = hi*exp(-0.255*frt*conc[npsi]);
	 hfo = ho*exp(0.255*frt*conc[npsi]);
       nadh= tnad - conc[nnad];
	 }
	 
void Ldistr::c3calc( double *py,double *pdydt) {
  double kr= nv.frw[rqp_FS]*hfo*hfo;
  double ox=0.955, bp=0.;
//Qo->FeS
double sum= pBC1q.shiftFeS(nv.frw[qp_FS],kr,nv.frw[qp_bl],nv.frw[rqp_bl], nv.frw[vros], ox, dconc[nc3ros]);
sum += qhpBC1.shiftFeS(nv.frw[qp_FS],kr,nv.frw[qp_bl],nv.frw[rqp_bl], nv.frw[vros], ox, dconc[nc3ros]);
bp = pBC1q.bypass(nv.frw[bypas]);
bp += qhpBC1.bypass(nv.frw[bypas]); sum -= bp;
//nv.flx[fros]=cros;
nv.flx[fbp]=bp;
		dconc[npsi] += sum*fc;
//1: e-trasport from FeS-protein (positions 3 or 1) to c1:
   sum=pBC1q.shift1(3,nv.frw[FS_c1],nv.frw[rFS_c1]);
   sum += qhpBC1.shift1(3,nv.frw[FS_c1],nv.frw[rFS_c1]);
   sum += BC1.shift1(1,nv.frw[FS_c1],nv.frw[rFS_c1]);
   sum += BC1qn.shift1(1,nv.frw[FS_c1],nv.frw[rFS_c1]);
		dconc[npsi] += sum*fc;

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
//		dconc[nhi] -= 2.*sum/buf/vol; 
  dconc[npsi] += 2.*sum*fc;
		
//6: binding qh2 on the p-side
  dconc[nqh] += BC1.qphbind(qhpBC1,conc[nqh],nv.frw[qHbnd],nv.frw[rqHbnd],bc15,1);
  dconc[nqh] += BC1qn.qphbind(pBC1q,conc[nqh],nv.frw[qHbnd],nv.frw[rqHbnd]); 

//9: dissociation of qh2 on the n-side
  dconc[nqh] += BC1qn.dissqn(BC1,conc[nqh],nv.frw[qhnds],nv.frw[rqhnds],bc15,1);
  dconc[nqh] += pBC1q.dissqn(qhpBC1,conc[nqh],nv.frw[qhnds],nv.frw[rqhnds]);

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
    dconc[npsi] += 2.*sum*fc;
    }
    
void Ldistr::tca(double *py,double *pdydt) {
//11: succinate dehydrogenase
double sum, vfac(1.);// = nv.flx[fsdh] = MM(nv.frw[vsdh],nv.nv[kq],qq)*MM(1.0,nv.nv[ksdh],conc[nsuc]);
//    dconc[nqh] += sum; 
//    dconc[nsuc] -= sum; 
//    dconc[nfum] += sum;
//Fumarate oxidation and MDH reaction
sum = nv.flx[fuoa] = nv.frw[vmdh]*(conc[nnad]*conc[nfum]-nadh*conc[noaa]);
  dconc[nnad] -= sum; 
  dconc[nfum] -= sum/vfac; 
  dconc[noaa] += sum/vfac; 
//pyruvate transport:
nv.flx[fpyrt] = nv.frw[vpyrTr]*(nv.nv[pyr_o]-conc[npyr]);
      dconc[npyr] += nv.flx[fpyrt];
// Citrate Synthase
sum =nv.flx[fcs] = nv.frw[vcs]*conc[noaa]*conc[npyr];
    dconc[noaa] -= sum/vfac;    
    dconc[npyr] -= sum/vfac;  
    dconc[ncit] += sum/vfac;  
      dconc[nnad] -= sum;

//TCA cycle from cit to suc
sum =nv.flx[ftca] = nv.frw[vtca]*conc[nnad]*conc[ncit];
    dconc[ncit] -= sum/vfac;   
    dconc[nnad] -= 2.*sum;  
    dconc[nsuc] += sum/vfac;  
// succinate exchange to fumarate/malate:
double a=0.001;
sum = nv.flx[fsm] =a*nv.frw[vsucmal]*(nv.nv[cmal_o]-conc[nfum]);
       dconc[nfum] += sum/vfac;
sum = nv.flx[fsm] =(1.-a)*nv.frw[vsucmal]*(nv.nv[suc_o]*conc[nfum]-conc[nsuc]*nv.nv[cmal_o]);
       dconc[nsuc] += sum/vfac;
       dconc[nfum] -= sum/vfac;
// succinate entry:
sum = nv.nv[vsuc_i]*(nv.nv[suc_o]-conc[nsuc]);
       dconc[nsuc] += sum/vfac;/**/
//malic enzyme:
sum = nv.flx[fme] = conc[nfum]*conc[nnad]*nv.frw[vme];
      dconc[nfum] -= sum/vfac;
      dconc[nnad] -= sum;
      dconc[npyr] += sum/vfac;
//leak
nv.flx[flk] = 0.005*conc[npsi]*leak(nv.frw[vlk]);
sum = nv.flx[flk];
    dconc[npsi] -= 2.*sum*fc;

   if(conc[nca]> 0.001) {
   a=conc[npsi]*frt;
   double e=exp(a);
   sum = 0.0005*a*conc[nca]*e/(1-e);
   dconc[nca] += sum;
   dconc[npsi] += 4.*sum*fc;
   }
    }

void Ldistr::c1calc( double *py,double *pdydt) {
//COMPLEX I (qp-qp-qn-qn-n2-fmn-fmn)
 nv.flx[fc1] = coreI.fmnred(nv.frw[vfred], nv.frw[vrfred], nadh,conc[nnad],fmn,fmnh);
 nv.flx[fc1] = cIq.fmnred(nv.frw[vfred], nv.frw[vrfred], nadh,conc[nnad],fmn,fmnh); 
 dconc[nnad] += nv.flx[fc1]; //NADH->FMN
 coreI.n562(nv.frw[vn56],nv.frw[vrn56],n5red, n6ar);
 cIq.n562(nv.frw[vn56],nv.frw[vrn56],n5red, n6ar);
	double kf1 = nv.frw[vn2qn1] * exp(-0.95*frt*conc[npsi]);
	double kr1 = nv.frw[vrn2qn1] * exp(0.95*frt*conc[npsi]); 
	double kf2 = nv.frw[vn2qn2];
	double kr2 = nv.frw[vrn2qn2]; 
 double sum = cIq.n2q(kf1,kr1,kf2,kr2,n2red);// 0.;// n2 -> Q
 dconc[npsi] += 4.*sum*fc; nv.flx[fc1] = sum;
 dconc[nqh] += cIq.qhdiss1(coreI,conc[nqh],nv.frw[vndis],nv.frw[vrndis]); // QH2->
 coreI.qbind1(cIq,qq,nv.frw[vpbind],nv.frw[vrpbind]);// QH2<-
}
void Ldistr::c2calc( double *py,double *pdydt) {
//COMPLEX I (qp-qp-qn-qn-n2-fmn-fmn)
 nv.flx[fc2] = coreII.fadred(nv.frw[vfadf], nv.frw[vfadr], conc[nsuc], conc[nfum]);
 nv.flx[fc2] = cIIq.fadred(nv.frw[vfadf], nv.frw[vfadr], conc[nsuc], conc[nfum]); 
 dconc[nfum] += nv.flx[fc2]/1000; 
 dconc[nsuc] -= nv.flx[fc2]/1000; 
	double kf1 = nv.frw[vbq1] ;
	double kr1 = nv.frw[vrbq1]; 
	double kf2 = nv.frw[vbq2];
	double kr2 = nv.frw[vrbq2];
 coreII.fadfs(kf1,kr1,f2s);
 cIIq.fadfs(kf1,kr1,f2s);
 double sum = cIIq.fsq(kf2,kr2,kf2,kr2,br);
 dconc[nqh] += cIIq.qhdiss1(coreII,conc[nqh],nv.frw[vqdis],nv.frw[vrqdis], cr11);
 coreII.qbind1(cIIq,qq,nv.frw[vqbind],nv.frw[vrqbind],cr11);
 dconc[nc2ros]=coreII.fadros(c2ros);
 dconc[nc2ros]+=cIIq.fadros(c2ros);
}

void Ldistr::distr( double *py,double *pdydt) {
  seteq(py,pdydt); 
  c3calc(py,pdydt);
  tca(py,pdydt);
  c1calc(py,pdydt);
  c2calc(py,pdydt);
}

