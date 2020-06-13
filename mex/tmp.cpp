//---------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include "nums.hh"
#include "modlab.h"
//---------------------------------------------------------------------------

using namespace std;
using namespace Label;
Ldistr horse;
//int  method=0;
int NN,numex; 
    const double tfac=5100.;

const int qp_FS(0), rqp_FS(qp_FS+1), FS_c1(rqp_FS+1), rFS_c1(FS_c1+1), qp_bl(rFS_c1+1), rqp_bl(qp_bl+1), bl_bh(rqp_bl+1), rbl_bh(bl_bh+1), bh_qn1(rbl_bh+1), rbh_qn1(bh_qn1+1), bh_qn2(rbh_qn1+1), rbh_qn2(bh_qn2+1), qHbnd(rbh_qn2+1), rqHbnd(qHbnd+1), qnbnd(rqHbnd+1), rqnbnd(qnbnd+1), qpdis(rqnbnd+1), rqpdis(qpdis+1), qhnds(rqpdis+1),rqhnds(qhnds+1), vc1c(rqhnds+1),bypas(vc1c+1), vsdh(bypas+1), vfadf(vsdh+1), vfadr(vfadf+1), vbq1(vfadr+1), vrbq1(vbq1+1), vbq2(vrbq1+1), vrbq2(vbq2+1), vqdis(vrbq2+1), vrqdis(vqdis+1), vqbind(vrqdis+1), vrqbind(vqbind+1), vfred(vrqbind+1), vrfred(vfred+1), vn56(vrfred+1), vrn56(vn56+1), vn2qn1(vrn56+1), vrn2qn1(vn2qn1+1), vn2qn2(vrn2qn1+1), vrn2qn2(vn2qn2+1), vndis(vrn2qn2+1), vrndis(vndis+1), vpbind(vrndis+1), vrpbind(vpbind+1), vros(vrpbind+1), vlk(vros+1), vatsyn(vlk+1), vtca(vatsyn+1), vmdh(vtca+1), vatase(vmdh+1), vox(vatase+1),vsuc_i(vox+1),vme(vsuc_i+1),vsuoa(vme+1),vpyrTr(vsuoa+1), vsucmal(vpyrTr+1), vcs(vsucmal+1), rct(vcs+1),ntmax(rct+1),ksdh(ntmax+1),kq(ksdh+1),cmal_o(kq+1), katase(cmal_o+1), pyr_o(katase+1), suc_o(pyr_o+1), oaa_o(suc_o+1), nNV(oaa_o+1);

//concentrations:
const int nqh(0), npsi(nqh+1), nnad(npsi+1), npyr(nnad+1), nsuc(npyr+1), nfum(nsuc+1), noaa(nfum+1), ncit(noaa+1), nc1c(ncit+1), nca(nc1c+1), numx(nca+1);
//fluxes:
const int fsdh(0), fsm(fsdh+1), fc1c(fsm+1), fme(fc1c+1), fpyrt(fme+1), fcs(fpyrt+1), ftca(fcs+1), fuoa(ftca+1), fatps(fuoa+1), flk(fatps+1), fbp(flk+1), fc1(fbp+1), fc2(fc1+1), flkc1(fc2+1), nflx(flkc1+1);

double Ldistr::fmnh[]={0, 9.056359999999999e-05, 0.000800642, 0.0065, 0.0989, 0.437, 0.804, 1};
double Ldistr::fs[]={0, 0.00146, 0.00433, 0.0122, 0.0447, 0.0712, 0.195, 0};
double Ldistr::fmn[]={1, 0.998, 0.994, 0.981, 0.856, 0.49, 0, 0};
double Ldistr::n5red[]={0, 0.249, 0.496, 0.739, 0.922, 0.97, 0.956, 1};
double Ldistr::n6ar[]={0, 0.0358, 0.501, 1};
double Ldistr::n2red[]={0, 0.928, 0.997, 1};
double Ldistr::f2s[]={0, 0.13511, 0.44567, 0.97849, 1};
double Ldistr::br[]={0, 0.29294, 0.68074, 0.99178, 1};
double Ldistr::bhr[]={0, 0.5486148322105801, 0.028990702336663, 1};
double Ldistr::sqi[]={0, 0.45138516778941, 0.028990702336663, 0};
double Ldistr::qhi[]={0, 0, 0.97100929766333, 1};
double Ldistr::fesr[]={0, 0.67727729241309, 1};
double Ldistr::c1r[]={0, 0.3227227075869, 1};

int Ldistr::setny() {
        pBC1q.ny = 0;
         pBC1q.ml = pBC1q.chkrl(6);
        BC1.ny = pBC1q.ml;
         BC1.ml = BC1.getlen();
          BC1.map=new int[BC1.ml];
          for(int i=0;i<BC1.ml;i++) BC1.map[i]=i;
        qhpBC1.ny = BC1.ny+BC1.ml; 
         qhpBC1.ml = qhpBC1.chk0();
        BC1qn.ny = qhpBC1.ny + qhpBC1.ml;
         BC1qn.ml = BC1qn.chkl(4);
        coreI.setny(BC1qn.ny + BC1qn.ml);   coreI.setnfmn2(31,1);
        cIq.setny(coreI.getny() + coreI.getlen()); cIq.setnfmn2(32,3);
        coreII.setnfsb(1); coreII.setny(cIq.ny + cIq.getlen());
        cIIq.setnfsb(3); cIIq.setny(coreII.getny() + coreII.getlen());
        nmet = cIIq.getny() + cIIq.getlen();   cout<<"nmet="<<nmet<<endl;/**/
        return nmet;
}
void Ldistr::setdisot(double *pyinit) {
        pBC1q.disot = &pyinit[pBC1q.ny];
        BC1.disot = &pyinit[BC1.ny];
        qhpBC1.disot = &pyinit[qhpBC1.ny];
        BC1qn.disot = &pyinit[BC1qn.ny];
        coreI.disot  = &pyinit[coreI.getny()];
        cIq.disot = &pyinit[cIq.getny()];
        coreII.disot  = &pyinit[coreII.getny()];
        cIIq.disot = &pyinit[cIIq.getny()];
        dconc = &pyinit[nmet];
}
void Ldistr::setisot(double *pyinit) {
        pBC1q.isot = &pyinit[pBC1q.ny];
        BC1.isot = &pyinit[BC1.ny];
        qhpBC1.isot = &pyinit[qhpBC1.ny];
        BC1qn.isot = &pyinit[BC1qn.ny];
        coreI.isot  = &pyinit[coreI.getny()];
        cIq.isot = &pyinit[cIq.getny()];
        coreII.isot  = &pyinit[coreII.getny()];
        cIIq.isot = &pyinit[cIIq.getny()];
        conc = &pyinit[nmet];
}

double nvv[]={8, 18, 350, 1500, 2500, 71, 900, 8.35, 5000, 7400, 5000000000000, 1000, 5.31969, 4.92438, 24.1563, 0.343086, 6.71485, 1.40863, 7.17725, 13.1662, 1.48978, 0, 2.58392, 0.05, 0, 25, 0, 35, 0, 2, 0, 3, 0, 1111, 4444, 8.93369, 15.2253, 116.351, 161.335, 973.389, 284.016, 41.4094, 12.7863, 21.625, 6.28562, 0, 0.05, 0, 1.6704, 0.963218, 0, 0.3, 0.143191, 6.33577e-06, 0, 0.070841, 0.00830269, 5.76002, 15, 4000, 0.5, 0.5, 0.3, 0, 0.005, 0.005, 0};
double ffrw[]={120, 270, 5250, 22500, 37500, 1065, 13500, 125.25, 75000, 111000, 75000000000000, 15000, 79.79535, 73.8657, 362.3445, 5.14629, 100.72275, 21.12945, 107.65875, 197.493, 22.3467, 0, 38.7588, 0.75, 0, 375, 0, 525, 0, 30, 0, 45, 0, 16665, 66660, 134.00535, 228.3795, 1745.265, 2420.025, 14600.835, 4260.240000000001, 621.141, 191.7945, 324.375, 94.2843, 0, 0.75, 0, 25.056, 14.44827, 0, 4.5, 2.147865, 9.503655e-05, 0, 1.062615, 0.12454035, 86.4003};

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
  double kr= ffrw[rqp_FS]*hfo*hfo;
  double cros=0., ox=0.955, bp=0.;
//Qo->FeS
double sum= pBC1q.shiftFeS(ffrw[qp_FS],kr,ffrw[qp_bl],ffrw[rqp_bl], ffrw[vros], ox, cros);
//sum += qhpBC1.shiftFeS(ffrw[qp_FS],kr,ffrw[qp_bl],ffrw[rqp_bl], ffrw[vros], ox, cros);
bp = pBC1q.bypass(ffrw[bypas]);
//bp += qhpBC1.bypass(ffrw[bypas]); sum -= bp;
//nv.flx[fros]=cros;
nv.flx[fbp]=bp;
		dconc[npsi] += sum*fc;
//1: e-trasport from FeS-protein (positions 3 or 1) to c1:
   sum=pBC1q.shift1(3,ffrw[FS_c1],ffrw[rFS_c1]);
//   sum += qhpBC1.shift1(3,ffrw[FS_c1],ffrw[rFS_c1]);
//   sum += BC1.shift1(1,ffrw[FS_c1],ffrw[rFS_c1]);
//   sum += BC1qn.shift1(1,ffrw[FS_c1],ffrw[rFS_c1]);
		dconc[npsi] += sum*fc;

//3: e-trasport between the hemes of cyt.b //bL->bH
double a(exp(-0.55*frt*conc[npsi])), b(exp(0.55*frt*conc[npsi]));
double kf= ffrw[bl_bh] * a;
       kr= ffrw[rbl_bh] * b;
              pBC1q.shift1(5,kf,kr);
//              qhpBC1.shift1(5,kf,kr);
//	      BC1.shift1(3,kf,kr);
//	      BC1qn.shift1(3,kf,kr);
	      
//4, 5: e-trasport further to q on n-side
	kf = ffrw[bh_qn2] * hfi * hfi;
  sum = pBC1q.bhqn(6,ffrw[bh_qn1],ffrw[rbh_qn1],kf,ffrw[rbh_qn2]);
//  sum += BC1qn.bhqn(4,ffrw[bh_qn1],ffrw[rbh_qn1],kf,ffrw[rbh_qn2]);
//		dconc[nhi] -= 2.*sum/buf/vol; 
		dconc[npsi] += 2.*sum*fc;
		
//6: binding qh2 on the p-side
  dconc[nqh] += BC1.qphbind(qhpBC1,conc[nqh],ffrw[qHbnd],ffrw[rqHbnd],bc15,1);
  dconc[nqh] += BC1qn.qphbind(pBC1q,conc[nqh],ffrw[qHbnd],ffrw[rqHbnd]); 

//9: dissociation of qh2 on the n-side
  dconc[nqh] += BC1qn.dissqn(BC1,conc[nqh],ffrw[qhnds],ffrw[rqhnds],bc15,1);
  dconc[nqh] += pBC1q.dissqn(qhpBC1,conc[nqh],ffrw[qhnds],ffrw[rqhnds]);

//7: binding of q on the n-side	
    BC1.bindqn(BC1qn,qq,ffrw[qnbnd],ffrw[rqnbnd],bc15,1);    
    qhpBC1.bindqn(pBC1q,qq,ffrw[qnbnd],ffrw[rqnbnd]);

//8: dissociation of q on the p-side 	
 qhpBC1.qpdiss(BC1,qq,ffrw[qpdis],ffrw[rqpdis]);
 pBC1q.qpdiss(BC1qn,qq,ffrw[qpdis],ffrw[rqpdis]);

//12: e-transport out from c1 (assumed to c)
        kf = ffrw[vc1c]*a;
	sum = pBC1q.c1c(4,kf, ox); 
//	sum += BC1.c1c(2,kf, ox,1);
//	sum += BC1qn.c1c(2,kf, ox);
//	sum += qhpBC1.c1c(4,kf, ox);
	nv.flx[fc1c] = sum;
    dconc[npsi] += 2.*sum*fc;
    }
    
void Ldistr::tca(double *py,double *pdydt) {
//11: succinate dehydrogenase
double sum;// = nv.flx[fsdh] = MM(ffrw[vsdh],nvv[kq],qq)*MM(1.0,nvv[ksdh],conc[nsuc]);
//    dconc[nqh] += sum; 
//    dconc[nsuc] -= sum; 
//    dconc[nfum] += sum;
//Fumarate oxidation and MDH reaction
sum = nv.flx[fuoa] = ffrw[vmdh]*(conc[nnad]*conc[nfum]-nadh*conc[noaa]);
  dconc[nnad] -= sum; 
  dconc[nfum] -= sum; 
  dconc[noaa] += sum; 
//pyruvate transport:
nv.flx[fpyrt] = ffrw[vpyrTr]*(nvv[pyr_o]-conc[npyr]);
      dconc[npyr] += nv.flx[fpyrt];
// Citrate Synthase
sum =nv.flx[fcs] = ffrw[vcs]*conc[noaa]*conc[npyr];
    dconc[noaa] -= sum;    
    dconc[npyr] -= sum;  
    dconc[ncit] += sum;  
      dconc[nnad] -= sum;

//TCA cycle from cit to suc
sum =nv.flx[ftca] = ffrw[vtca]*conc[nnad]*conc[ncit];
    dconc[ncit] -= sum;   
    dconc[nnad] -= 2.*sum;  
    dconc[nsuc] += sum;  
// succinate exchange to fumarate/malate:
double a=0.001;
sum = nv.flx[fsm] =a*ffrw[vsucmal]*(nvv[cmal_o]-conc[nfum]);
       dconc[nfum] += sum;
sum = nv.flx[fsm] =(1.-a)*ffrw[vsucmal]*(nvv[suc_o]*conc[nfum]-conc[nsuc]*nvv[cmal_o]);
       dconc[nsuc] += sum;
       dconc[nfum] -= sum;
// succinate entry:
sum = nvv[vsuc_i]*(nvv[suc_o]-conc[nsuc]);
       dconc[nsuc] += sum;/**/
//malic enzyme:
sum = nv.flx[fme] = conc[nfum]*conc[nnad]*ffrw[vme];
      dconc[nfum] -= sum;
      dconc[nnad] -= sum;
      dconc[npyr] += sum;
//leak
nv.flx[flk] = 0.005*conc[npsi]*leak(ffrw[vlk]);
sum = nv.flx[flk];
    dconc[npsi] -= 2.*sum*fc;

   if(conc[nca]> 0.1) {
   a=dconc[npsi]*frt;
   double e=exp(a);
   sum = 0.0025*a*conc[nca]*e/(1-e);
   dconc[nca] += sum;
   dconc[npsi] += 4.*sum*fc;
   }
    }

void Ldistr::c1calc( double *py,double *pdydt) {
//COMPLEX I (qp-qp-qn-qn-n2-fmn-fmn)
//nv.flx[fc1] = coreI.fmnred(ffrw[vfred], ffrw[vrfred], nadh,conc[nnad],fmn,fmnh);
nv.flx[fc1] = cIq.fmnred(ffrw[vfred], ffrw[vrfred], nadh,conc[nnad],fmn,fmnh); 
dconc[nnad] += nv.flx[fc1]; //NADH->FMN
// coreI.n562(ffrw[vn56],ffrw[vrn56],n5red, n6ar);
 cIq.n562(ffrw[vn56],ffrw[vrn56],n5red, n6ar);
	double kf1 = ffrw[vn2qn1] * exp(-0.95*frt*conc[npsi]);
	double kr1 = ffrw[vrn2qn1] * exp(0.95*frt*conc[npsi]); 
	double kf2 = ffrw[vn2qn2];
	double kr2 = ffrw[vrn2qn2]; 
  double sum = cIq.n2q(kf1,kr1,kf2,kr2,n2red);// 0.;// n2 -> Q
dconc[npsi] += 4.*sum*fc; nv.flx[fc1] = sum;
dconc[nqh] += cIq.qhdiss1(coreI,conc[nqh],ffrw[vndis],ffrw[vrndis]); // QH2->
 coreI.qbind1(cIq,qq,ffrw[vpbind],ffrw[vrpbind]);// QH2<-
}
void Ldistr::c2calc( double *py,double *pdydt) {
//COMPLEX I (qp-qp-qn-qn-n2-fmn-fmn)
// nv.flx[fc2] = coreII.fadred(ffrw[vfadf], ffrw[vfadr], conc[nsuc], conc[nfum]);
 nv.flx[fc2] = cIIq.fadred(ffrw[vfadf], ffrw[vfadr], conc[nsuc], conc[nfum]); 
 dconc[nfum] += nv.flx[fc2]; 
 dconc[nsuc] -= nv.flx[fc2]; 
	double kf1 = ffrw[vbq1] ;
	double kr1 = ffrw[vrbq1]; 
	double kf2 = ffrw[vbq2];
	double kr2 = ffrw[vrbq2];
// coreII.fadfs(kf1,kr1,f2s);
 cIIq.fadfs(kf1,kr1,f2s);
 double sum = cIIq.fsq(kf2,kr2,kf2,kr2,br);
 dconc[nqh] += cIIq.qhdiss1(coreII,conc[nqh],ffrw[vqdis],ffrw[vrqdis], cr11);
 coreII.qbind1(cIIq,qq,ffrw[vqbind],ffrw[vrqbind],cr11);
// coreII.fadros(0.00001);
// cIIq.fadros(0.00001);
}

void Ldistr::distr( double *py,double *pdydt,double p1) {
  nvv[suc_o]=p1;
  seteq(py,pdydt); 
  c3calc(py,pdydt);
  tca(py,pdydt);
  c1calc(py,pdydt);
  c2calc(py,pdydt);
  delete[] BC1.map;
  delete[] qhpBC1.map;
  delete[] BC1qn.map;
  delete[] pBC1q.map;
}

int main(){
  int nmet=horse.setny(); NN=numx+nmet;
  double yprime[NN], py[NN];
  ifstream fn("i2");
  for(int i=0;i<NN;i++) fn>>py[i];
  horse.distr(py, yprime,0.005);
  for(int i=0;i<NN;i++) cout<<i<<") "<<py[i]<<"  "<<yprime[i]<<'\n';
  return 0;
}

