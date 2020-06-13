//---------------------------------------------------------------------------
#include "mex.h"
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
    const double tfac=2100.;

const int qp_FS(0), rqp_FS(qp_FS+1), FS_c1(rqp_FS+1), rFS_c1(FS_c1+1), qp_bl(rFS_c1+1), rqp_bl(qp_bl+1), bl_bh(rqp_bl+1), rbl_bh(bl_bh+1), bh_qn1(rbl_bh+1), rbh_qn1(bh_qn1+1), bh_qn2(rbh_qn1+1), rbh_qn2(bh_qn2+1), qHbnd(rbh_qn2+1), rqHbnd(qHbnd+1), qnbnd(rqHbnd+1), rqnbnd(qnbnd+1), qpdis(rqnbnd+1), rqpdis(qpdis+1), qhnds(rqpdis+1),rqhnds(qhnds+1), vc1c(rqhnds+1),bypas(vc1c+1), vsdh(bypas+1), vfred(vsdh+1), vrfred(vfred+1), vfn2(vrfred+1), vrfn2(vfn2+1), vn2qn1(vrfn2+1), vrn2qn1(vn2qn1+1), vqpqn(vrn2qn1+1), vrqpqn(vqpqn+1), vndis(vrqpqn+1), vrndis(vndis+1), vpbind(vrndis+1), vrpbind(vpbind+1), vn2qn2(vrpbind+1), vrn2qn2(vn2qn2+1), vros(vrn2qn2+1), vlk(vros+1), vatsyn(vlk+1), vtca(vatsyn+1), vmdh(vtca+1), vatase(vmdh+1), vox(vatase+1),vsuc_i(vox+1),vme(vsuc_i+1),vsuoa(vme+1),vpyrTr(vsuoa+1), vsucmal(vpyrTr+1), vcs(vsucmal+1), rct(vcs+1),ntmax(rct+1),ksdh(ntmax+1),kq(ksdh+1),cmal_o(kq+1), katase(cmal_o+1), pyr_o(katase+1), suc_o(pyr_o+1), oaa_o(suc_o+1), nNV(oaa_o+1);
//concentrations:
const int nqh(0), npsi(nqh+1), nnad(npsi+1), npyr(nnad+1), nsuc(npyr+1), nfum(nsuc+1), noaa(nfum+1), ncit(noaa+1), nc1c(ncit+1), numx(nc1c+1);
//fluxes:
const int fsdh(0), fsm(fsdh+1), fsuoa(fsm+1), fme(fsuoa+1), fpyrt(fme+1), fcs(fpyrt+1), ftca(fcs+1), fuoa(ftca+1), fatps(fuoa+1), flk(fatps+1), fatpa(flk+1), fc1(fatpa+1), flkc1(fc1+1), nflx(flkc1+1);

int Ldistr::setny() {
        pBC1q.ny = 0;
         pBC1q.ml = pBC1q.chkrl(6);
        BC1.ny = pBC1q.ml;
         BC1.ml = 16;
          BC1.map=new int[BC1.ml];
          for(int i=0;i<BC1.ml;i++) BC1.map[i]=i;
        qhpBC1.ny = BC1.ny+BC1.ml; 
         qhpBC1.ml = qhpBC1.chk0();
        BC1qn.ny = qhpBC1.ny + qhpBC1.ml;
         BC1qn.ml = BC1qn.chkl(4);
        qn1.ny   = BC1qn.ny + BC1qn.ml;
        qnp1.ny  = qn1.ny + 1;
        nmet = qnp1.ny + qnp1.getlen(); 
        return nmet;
}
void Ldistr::setdisot(double *pyinit) {
        pBC1q.disot = &pyinit[pBC1q.ny];
        BC1.disot = &pyinit[BC1.ny];
        qhpBC1.disot = &pyinit[qhpBC1.ny];
        BC1qn.disot = &pyinit[BC1qn.ny];
        qn1.disot  = &pyinit[qn1.ny];
        qnp1.disot = &pyinit[qnp1.ny];
        dconc = &pyinit[nmet];
}
void Ldistr::setisot(double *pyinit) {
        pBC1q.isot = &pyinit[pBC1q.ny];
        BC1.isot = &pyinit[BC1.ny];
        qhpBC1.isot = &pyinit[qhpBC1.ny];
        BC1qn.isot = &pyinit[BC1qn.ny];
        qn1.isot  = &pyinit[qn1.ny];
        qnp1.isot = &pyinit[qnp1.ny];
        conc = &pyinit[nmet];
}


void Ldistr::distr( double *py,double *pdydt,double p1,double p2) {
double nvv[]={8.614000000000001, 83.5274, 271.783, 1261.108, 1310.618, 0.0035, 400.024, 2.20189, 172.421, 257.1712, 56788300000, 251.6056, 5.31969, 4.92438, 24.1563, 0.343086, 6.71485, 1.40863, 7.17725, 13.1662, 1.489785, 0.0005, 2, 702, 1500, 150, 0.06, 950, 16, 6029820000, 22003500, 406.946, 38.8533, 382.912, 5.38921, 30000000000, 9, 0.00010475, 2.2, 0, 1.66841, 0.861088, 0, 0.217603, 0.715709, 4.37521e-07, 0.00084913, 0.1, 0.013, 4.01673, 15, 4000, 0.5, 0.5, 0, 0, 0.02, 0.0255, 0};
double ffrw[]={129.21, 1252.911, 4076.745, 18916.62, 19659.27, 0.0525, 6000.36, 33.02835, 2586.315, 3857.568, 851824500000, 3774.084, 79.79535, 73.8657, 362.3445, 5.14629, 100.72275, 21.12945, 107.65875, 197.493, 22.346775, 0.0075, 30, 10530, 22500, 2250, 0.8999999999999999, 14250, 240, 90447300000, 330052500, 6104.190000000001, 582.7995, 5743.679999999999, 80.83815, 450000000000, 135, 0.00157125, 33, 0, 25.02615, 12.91632, 0, 3.264045, 10.735635, 6.562815e-06, 0.01273695, 1.5, 0.195, 60.25095};
nvv[suc_o]=p1;
int nmet=setny(); NN=numx+nmet;
	setisot(py);
	setdisot(pdydt);
    pBC1q.total(); qhpBC1.total(); BC1qn.total(); BC1.total();
    qnp1.total(); 
 double ciii0=c3t-pBC1q.cont-qhpBC1.cont-BC1qn.cont-BC1.cont;
 double cin11000 = c1t-qnp1.cont-qn1.isot[in10100];
 double qq=qt-2*(pBC1q.cont+ciii0)-qhpBC1.cont-BC1qn.cont-conc[nqh]-cin11000-qn1.isot[0]-2.*qnp1.cont;
	for (int i=0;i<NN;i++) pdydt[i]=0.;
	double hfi = hi*exp(-0.255*frt*conc[npsi]);
	double hfo = ho*exp(0.255*frt*conc[npsi]);
//COMPLEX III
//0, 2, 10: qh2-bound -> FeS (2->3) and then to cytB (1->5), ROS production:
//proton transport & potential change
  double kr= ffrw[rqp_FS]*hfo*hfo;
  double cros=0., ox=0.955, bp=0.;
double sum= pBC1q.shiftFeS(ffrw[qp_FS],kr,ffrw[qp_bl],ffrw[rqp_bl],ffrw[vros],ox,cros);
sum += qhpBC1.shiftFeS(ffrw[qp_FS],kr,ffrw[qp_bl],ffrw[rqp_bl],ffrw[vros],ox, cros);
nv.flx[fatps]=cros;
bp = pBC1q.bypass(ffrw[bypas]);
bp += qhpBC1.bypass(ffrw[bypas]); sum -= bp;
nv.flx[fatpa]=bp;
//		dconc[nho] += 2.*sum/buf/vol; 
		dconc[npsi] += sum*fc;
//1: e-trasport from FeS-protein (positions 3 or 1) to c1:
   sum=pBC1q.shift1(3,ffrw[FS_c1],ffrw[rFS_c1]);
   sum += qhpBC1.shift1(3,ffrw[FS_c1],ffrw[rFS_c1]);
   sum += BC1.shift1(1,ffrw[FS_c1],ffrw[rFS_c1]);
   sum += BC1qn.shift1(1,ffrw[FS_c1],ffrw[rFS_c1]);
		dconc[npsi] += sum*fc;
//3: e-trasport between the hemes of cyt.b   
double a(exp(-0.55*frt*conc[npsi])), b(exp(0.55*frt*conc[npsi]));
double kf= ffrw[bl_bh] * a;
       kr= ffrw[rbl_bh] * b;
              pBC1q.shift1(5,kf,kr);
              qhpBC1.shift1(5,kf,kr);
	      BC1.shift1(3,kf,kr);
	      BC1qn.shift1(3,kf,kr);
//4, 5: e-trasport further to q on n-side
	kf = ffrw[bh_qn2] * hfi * hfi;
  sum = pBC1q.bhqn(6,ffrw[bh_qn1],ffrw[rbh_qn1],kf,ffrw[rbh_qn2]);
  sum += BC1qn.bhqn(4,ffrw[bh_qn1],ffrw[rbh_qn1],kf,ffrw[rbh_qn2]);
//		dconc[nhi] -= 2.*sum/buf/vol; 
		dconc[npsi] += 2.*sum*fc;
//6: binding qh2 on the p-side
       dconc[nqh] += BC1.qphbind(qhpBC1,conc[nqh],ffrw[qHbnd],ffrw[rqHbnd]);
       dconc[nqh] += BC1qn.qphbind(pBC1q,conc[nqh],ffrw[qHbnd],ffrw[rqHbnd]); 
//9: dissociation of qh2 on the n-side
	dconc[nqh] += BC1qn.dissqn(BC1,conc[nqh],ffrw[qhnds],ffrw[rqhnds]);
	dconc[nqh] += pBC1q.dissqn(qhpBC1,conc[nqh],ffrw[qhnds],ffrw[rqhnds]);
//7: binding of q on the n-side	
    BC1.bindqn(BC1qn,qq,ciii0,ffrw[qnbnd],ffrw[rqnbnd],0);    
    qhpBC1.bindqn(pBC1q,qq,ciii0,ffrw[qnbnd],ffrw[rqnbnd],1);
//8: dissociation of q on the p-side 	
 qhpBC1.qpdiss(BC1,qq,ciii0,ffrw[qpdis],ffrw[rqpdis],0);
 pBC1q.qpdiss(BC1qn,qq,ciii0,ffrw[qpdis],ffrw[rqpdis],1);
//12: e-transport out from c1 (assumed to c)
        kf = ffrw[vc1c]*a;
	sum = pBC1q.c1c(4,kf, ox); 
	sum += BC1.c1c(2,kf, ox);
	sum += BC1qn.c1c(2,kf, ox);
	sum += qhpBC1.c1c(4,kf, ox);
	nv.flx[fsuoa] = sum;
//    dconc[nhi] -= 4.*sum/buf/vol; 
//    dconc[nho] += 2.*sum/buf/vol; 
    dconc[npsi] += 4.*sum*fc;
    dconc[nc1c] = nv.flx[fsuoa]-conc[nc1c];
//    dconc[nprod] += sum;   dconc[nprod] =0.; 
 // dconc[nox] += ffrw[vox]*(1.0-conc[nox]) - sum;
 //     dconc[nox] += 2.*o2deliv(conc[nox])/tfac - sum;
    
//11: succinate dehydrogenase
sum = nv.flx[fsdh] = MM(ffrw[vsdh],nvv[kq],qq)*MM(1.0,nvv[ksdh],conc[nsuc]);
    dconc[nqh] += sum; 
    dconc[nsuc] -= sum; 
    dconc[nfum] += sum;
double nadh= tnad - conc[nnad];
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
sum = nv.flx[fsm] =ffrw[vsucmal]*(nvv[suc_o]*conc[nfum]-conc[nsuc]*nvv[cmal_o]);
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
nv.flx[flk] = leak(ffrw[vlk]);
sum = nv.flx[flk];
//    dconc[nhi] += sum/buf/vol; 
//    dconc[nho] -= sum/buf/vol; 
    dconc[npsi] -= 2.*sum*fc;
//    dconc[ncit] = 0.0;
//    dconc[noaa] = 0.0;
    
//COMPLEX I (qp-qp-qn-qn-n2-fmn-fmn)
nv.flx[fc1] = qnp1.fmnred(hi,i1100000,i1100011, ffrw[vfred], ffrw[vrfred], nadh,conc[nnad]); //0:NADH->FMN
dconc[nnad] += nv.flx[fc1]; //0:NADH->FMN
  qnp1.redox1(i1100011,i1100101,ffrw[vfn2],ffrw[vrfn2]);     //1: FMN -> N2 (first e)
  qnp1.redox1(i1100101,i1101001,ffrw[vn2qn1],ffrw[vrn2qn1]); //2: N2 -> Qn (first e) 
	kf = ffrw[vqpqn] * hfi * hfi;
	kr = ffrw[vrqpqn] * hfo * hfo; 
  sum = qnp1.redox1(i1101001,i1011001,kf,kr);//3: interaction Qp & Qn
  qnp1.redox1(i1011001,i1110001,ffrw[vndis],ffrw[vrndis]);//4: exchange Qp <-> Qn
qnp1.redox1(i1110001,i1110100,ffrw[vfn2],ffrw[vrfn2]);    //6: FMN -> N2 (second e)
  sum += qnp1.redox1(i1110100,i1011100,kf,kr);//7: Qp -> Qn 
  nv.flx[flkc1]=sum;
//dconc[nho] += 2.*sum/buf/vol;
//dconc[nho] = 0.;
//dconc[nhi] -= 2.*sum/buf/vol; 
dconc[npsi] += 4.*sum*fc;
dconc[nqh] += qnp1.qpdiss1(qn1,conc[nqh],i1011100,in10100,ffrw[vndis],ffrw[vrndis]); //8: QH2->
	kf = ffrw[vn2qn2]*hfi*hfi;
  sum = qn1.redox2(in10100,cin11000,kf,ffrw[vrn2qn2]);//9: N2 -> Qn (second e)
//dconc[nhi] -= 2.*sum/buf/vol; 
//dconc[nhi] =0.; 
dconc[npsi] += 2.*sum*fc;
 qn1.qpbind1(qnp1,qq,cin11000,i1100000,ffrw[vpbind],ffrw[vrpbind]);//12: QH2<-
//	 	      if ((qq < -0.5)||(qq > 5)) throw "--negative conc--";
//  ROS in complex I: i1100000,i1100101, i1100011, i1110100, i1110001, i1101001, i1011100, i1011001;
   kf= ffrw[vros]*ox;
  sum = qnp1.redox1(i1101001,i1100000,kf,0.);
  sum += qnp1.redox1(i1110100,i1100000,kf,0.);
  sum += qnp1.redox1(i1110001,i1100000,kf,0.);
  sum += qnp1.redox1(i1101001,i1100000,kf,0.);
/**/delete[] pBC1q.map;
delete[] BC1.map;
delete[] qhpBC1.map;
delete[] BC1qn.map;
}

void mexFunction (int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
  int ip1, ip2; 
  double *par1,*par2, *py, *pdydt;

       plhs[0] = (mxArray *) mxCreateNumericArray 
        (2, mxGetDimensions(prhs[0]), mxDOUBLE_CLASS,mxREAL);
       py = mxGetPr (prhs[0]);
       pdydt = mxGetPr (plhs[0]);
       par1= mxGetPr(prhs[1]);
       par2= mxGetPr(prhs[2]);/**/
     horse.distr(py,pdydt,par1[0],par2[0]);
      }

