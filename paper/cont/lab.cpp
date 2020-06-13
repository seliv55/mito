//---------------------------------------------------------------------------
#include <fstream>
#include <iostream>
#include <iomanip>
#include "nr.h"
#include "nums.hh"
#include "modlab.h"
#include "main.hh"
//---------------------------------------------------------------------------
using namespace std;
using namespace Label;

void Ldistr::read (string fn) {
   char aa[20];
    ifstream fi(fn.c_str());
      pBC1q.read(fi); // pBC1q.lmpl(6);//fi >>aa;
      BC1.read(fi);
      qhpBC1.read(fi);//fi <<endl<<"qhpBC1/ ";  
      BC1qn.read(fi);
     coreI.read(fi);//fi <<endl<<"cIqnp/ ";   
     cIq.read(fi);//fi <<endl<<"cIqnp/ ";   
     coreII.read(fi);//fi <<endl<<"cIqnp/ ";   
     cIIq.read(fi);//fi <<endl<<"cIqnp/ ";   
//     c3t=pBC1q.cont+BC1.cont+qhpBC1.cont+BC1qn.cont+bc15;
      cout<<"c3t="<<c3t<<endl;
	fi >>aa>> conc[nqh];
	fi >>aa>> conc[npsi];
	fi >>aa>> conc[nnad];
	fi >>aa>> conc[npyr];
	fi >>aa>> conc[nsuc];
	fi >>aa>> conc[nfum];
	fi >>aa>> conc[noaa];
	fi >>aa>> conc[ncit];
	fi >>aa>> conc[nc1c];
	conc[nca]=0.; conc[nc2ros]=0.; conc[nc3ros]=0.;
	cout<<"psi="<<conc[npsi]<<'\n';
}
int Ldistr::readexp (int col, string fn) {
   ifstream fi(fn.c_str()); int i=0,j; string aaa;
        getline(fi,aaa); 
   for(j=0;j<col;j++) fi>>tshift; getline(fi,aaa); 
   while(!fi.eof()) { 
     fi>>tex[i]; tex[i] -= tshift; 
      if(tex[i]>=0.){for(j=0;j<col;j++) fi>>ex1[j][i]; i++;}
        getline(fi,aaa); } 
      fi.close(); return iter= i - 2;}

void Ldistr::write (string fn) const {
  ofstream fi(fn.c_str());
   fi.precision(16);
      pBC1q.write(fi);fi <<"\n"; //fi <<endl<<"pBC1q/ ";  
      BC1.write(fi);fi <<"\n"; 
      qhpBC1.write(fi);fi <<"\n"; //fi <<endl<<"qhpBC1/ ";  
      BC1qn.write(fi);fi <<"\n"; //fi <<endl<<"BC1qn/ ";  
     coreI.write(fi);fi <<"\n"; 
     cIq.write(fi);fi <<"\n"; //fi <<endl<<"cIqnp/ ";   
     coreII.write(fi);fi <<"\n";//fi <<endl<<"cIqnp/ ";   
     cIIq.write(fi);fi <<"\n";//fi <<endl<<"cIqnp/ ";   

	      fi <<"qh= "<< conc[nqh]<<"\n";// 
	      fi <<"psi= "<< conc[npsi]<<"\n";// 
	      fi <<"nad= "<< conc[nnad]<<"\n";// 
	      fi <<"pyr= "<< conc[npyr]<<"\n";// 
	      fi <<"suc= "<< conc[nsuc]<<"\n";// 
	      fi <<"fum= "<< conc[nfum]<<"\n";// 
	      fi <<"oa= "<< conc[noaa]<<"\n";// 
	      fi <<"cit= "<< conc[ncit]<<"\n";// 
	      fi <<"c1c= " << conc[nc1c]<<"\n";//flx 
//	      fi <<" "<< conc[nadp];//adp= 
}

void Ldistr::setarrays() {
  string sss;
  ifstream fi("c1rocalc");
  if(fi.is_open()){
    getline(fi,sss);
    for(int i=1;i<7;i++) 
       fi>>sss>>fmnh[i]>>fs[i]>>fmn[i]>>sss>>sss >> sss >> n5red[i];
    fi.close(); }
  else cout<<"unable to open \"c1rocalc\""<<endl;
  fmn[0]=1.; fmnh[0]=fs[0]=n5red[0]=0.;
  fmn[7]=fs[7]=0.; fmnh[7]=n5red[7]=1.;
  n6ar[0]=n2red[0]=0.;
  n6ar[1]=n2red[1]=1.;

  fi.open("c2rocalc");
  if(fi.is_open()){
    getline(fi,sss); 
    for(int i=0;i<5;i++) fi>>f2s[i]>>sss>>sss>>br[i];
    fi.close();}
  else cout<<"unable to open \"c2rocalc\""<<endl;

  fi.open("c3redox");
  if(fi.is_open()){
    getline(fi,sss); 
    for(int i=0;i<4;i++) fi>>bhr[i]>>sqi[i]>>qhi[i];
    getline(fi,sss); getline(fi,sss); 
    for(int i=0;i<3;i++) fi>>fesr[i]>>c1r[i];
    fi.close(); }
  else cout<<"unable to open \"c3redox\""<<endl;
//  ofstream fo("arrays");
//  fo.precision(16);
//complex1
//    fo<<"horse.fmnh[]={"<<fmnh[0];
//    for(int i=1;i<=7;i++) fo<<", "<<fmnh[i]; fo<<"};\n";
//    fo<<"horse.fs[]={"<<fs[0];
//    for(int i=1;i<=7;i++) fo<<", "<<fs[i]; fo<<"};\n";
//    fo<<"horse.fmn[]={"<<fmn[0];
//    for(int i=1;i<=7;i++) fo<<", "<<fmn[i]; fo<<"};\n";
//    fo<<"horse.n5red[]={"<<n5red[0];
//    for(int i=1;i<=7;i++) fo<<", "<<n5red[i]; fo<<"};\n";
//    fo<<"horse.n6ar[]={"<<n6ar[0];
//    for(int i=1;i<=3;i++) fo<<", "<<n6ar[i]; fo<<"};\n";
//    fo<<"horse.n2red[]={"<<n2red[0];
//    for(int i=1;i<=3;i++) fo<<", "<<n2red[i]; fo<<"};\n";
////complex2
//    fo<<"horse.f2s[]={"<<f2s[0];
//    for(int i=1;i<5;i++) fo<<", "<<f2s[i]; fo<<"};\n";
//    fo<<"horse.br[]={"<<br[0];
//    for(int i=1;i<5;i++) fo<<", "<<br[i]; fo<<"};\n";
////complex3
//    fo<<"horse.bhr[]={"<<bhr[0];
//    for(int i=1;i<4;i++) fo<<", "<<bhr[i]; fo<<"};\n";
//    fo<<"horse.sqi[]={"<<sqi[0];
//    for(int i=1;i<4;i++) fo<<", "<<sqi[i]; fo<<"};\n";
//    fo<<"horse.qhi[]={"<<qhi[0];
//    for(int i=1;i<4;i++) fo<<", "<<qhi[i]; fo<<"};\n";
//    fo<<"horse.fesr[]={"<<fesr[0];
//    for(int i=1;i<3;i++) fo<<", "<<fesr[i]; fo<<"};\n";
//    fo<<"horse.c1r[]={"<<c1r[0];
//    for(int i=1;i<3;i++) fo<<", "<<c1r[i]; fo<<"};\n";

//    fi.close();
} 
 
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
        coreI.setny(BC1qn.ny + BC1qn.ml);   coreI.setnfmn2(15,1);
        cIq.setny(coreI.getny() + coreI.getlen()); cIq.setnfmn2(16,3);
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

double Ldistr::fiout(double t, ostringstream& fi,int ii){
 double tmod = t/tfac;
sp3 = pBC1q.percent(3,1)+qhpBC1.percent(3,1);// 
   fsq = coreI.getfs(fs) + cIq.getfs(fs);
   sq1 = cIq.getsq();
   double fsc2= coreII.getfs() + cIIq.getfs(); 
   double qsc2= cIIq.getsq();
//   ros=sp3+fsq+sq1;
         double ai=3.-log10(hi);
if(!ii) fi<<"Time(min)"<<" prod"<<" sp3"<<" fsc2" <<" qsc2" <<" fum"<<" Bypass"<< " qh" << " psi" << " fc1" <<" ROS_c2" << " ROS_c3"<<" fc2"<<endl;
fi <<tmod<<" " << nv.flx[fc1c]*tfac/60./0.4 <<" "<< sp3/c3t <<" "<< fsc2/c2t <<" " << qsc2 <<" " << conc[nfum] <<" "<< nv.flx[fbp]*tfac/60./0.4 <<" " <<  (conc[nqh]) << " " << conc[npsi] << " " << nv.flx[fc1]*tfac/60./0.4 << " " << conc[nc2ros]/c2t <<" "<< conc[nc3ros]/c3t <<" "<<nv.flx[fc2]*tfac/60./0.4<<'\n';
     return conc[npsi];
}

void Ldistr::tout(int ipar,int& ii, ofstream& fi){
  static double t0=0., a=0.;
   if(!ii) {fi<<"par"<<" time "<< xi1 <<" "<< xi2 <<endl; ii++;}
   else fi << ipar <<" "<< tf <<" "<< xi1 <<" " << xi2 <<" "<< nv.nv[ipar] <<endl;
}/**/

void Ldistr::priout(){
static bool k=0;
if(!k) cout <<setw(9)<<"sq@Qo"<<setw(9)<<"nad"<<setw(9)<<"suc"<<setw(9)<<"sq@c1"<<setw(9)<<"prod"<< setw(9)<< "fmn" << setw(7) <<"psi" << setw(9) << "oaa" << setw(9) << "qh" <<  endl;

cout <<setw(9)<< sp3 <<setw(9)<< 0.<<setw(9)<< 0. <<setw(9)<< sn1 <<setw(9)<<dpdt <<setw(9)<< fmn << setw(7) << conc[npsi] << setw(9) << 0. <<setw(9)<< 0. <<setw(9)<< conc[nqh]<<endl;
k=1;
}

