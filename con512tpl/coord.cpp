#include <iostream>
#include <iomanip>
#include <sstream>
#include "nr.h"
#include "nums.hh"
#include "modlab.h"
#include "main.hh"
#include "analis.h"
using namespace std;
void Analis::descent(double *py,double nv2[],int par[], int npf){ 
   const double xili(0.9999);
   int ii(0),jfail(0), ipar;
     double factor(1.1), b, a, aa;
	while (npf>0)  {cout<<"npf="<<npf;
	ipar = rand() % npf; npf--; 
	int flag = 0; cout<<"; ipar="<<ipar<<"; par="<<par[ipar]<<endl;
   while (flag<2) { 
     a=horse.nv.changeval(par[ipar],factor);
	try{ 
xi=horse.efit(py);
b = (double)tf/(double)tfirst;
} catch( char const* str ) { cout << "exception: "<< str <<endl; xi=xi0; b=11;	}
 if (((xi*b)<xi0)&&(ii<3))   {
          horse.nv.nv2st(nv2);//saves nv&xx
          xi0 = xi; tfirst=tf; ii++;
                        horse.nv.write(ifn,xi,1);
   cout<<par[ipar]<<": xi="<<xi<<"; b="<<b<<endl;
       	}
                else {factor = 1./factor; flag++; ii=0;
	horse.nv.setval(par[ipar],a);//gets nv&xx
                 }
        }//end while flag
	if(npf) {for (int k=ipar;k<npf;k++) par[k] = par[k+1]; cout<<" par[npf]="<<par[npf-1]<<endl;}
	}
}

void Analis::coord(double ystart[]){
  const double f1(0.25);   double fact, nv1[nNV],nv2[nNV];
       horse.nv.nv2st(nv1);//saves nv&xx
cout<<"\nPerturbation+Coordinate Descent:\nPar#   xi2con  xi2iso  ParValue:"<<endl;
for(;;){ if (ifn>9000) return;	int sign,i=0;
 int par[]={22, 25, 26, 27, 28, 29, 30, 31, 32, 34, 36, 38, 39, 42, 43, 45, 46, 47, -18};
 //1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18, 19,20,21,
	while(par[i]>=0){
	 sign = rand() % 100;
	   fact = 1.- f1*(0.5 - sign*0.01); 
	  horse.nv.changeval(par[i], fact); i++;
	} cout<<par[i]<<"; last par; "<<endl;
	double *py = &ystart[0];
	try {
		xi0=horse.efit(py); tfirst=tf; xi=0.; cout<<"xi0="<<xi0<<endl;
while (xi < xi0) { 
  horse.nv.nv2st(nv2);//saves nv&xx
   descent(py,nv2,par,i);
      horse.nv.write(ifn,xi,1);
                xi=horse.efit(py);
      }
 } catch( char const* str ){cout << "exception: "<< str <<endl;	horse.nv.st2nv(nv1);}
	}
}
/*
void Analis::sensitivity(Vec_DP &ystart){
   const double f1(1.01);
   double xi01=x01, xi02=x02, xi03=x03,xi04=x04, nv1[nNV],st40[4],rot0[4],st30[4],st3pyr0[4],st4pyr0[4];
   for(int j=0;j<4;j++) {st40[j]=st4[j]; rot0[j]=rot[j]; st30[j]=st3[j]; st3pyr0[j]=st3pyr[j]; st4pyr0[j]=st4pyr[j];}
   int i=0,*ip=&horse.nv.par[1];
   ofstream fi("sensitivity.csv");
     cout<<"\nsensitivity"<<endl;
       horse.nv.nv2st(nv1);//saves nv&xx
	while(ip[i]>=0){
	 double a=horse.nv.changeval(ip[i], f1);
	double *py = &ystart[0];
  try {	comb(py);
cout<<ip[i]<<") "<<horse.nv.pname[ip[i]]<<": "<<(x01-xi01)<<", "<<(x02-xi02)<<", "<<(x03-xi03)<<", "<<(x04-xi04)<<endl;
fi<<ip[i]<<") "<<horse.nv.pname[ip[i]]<<": "<<(x01-xi01)<<", "<<(x02-xi02)<<", "<<(x03-xi03)<<", "<<(x04-xi04);
fi<<" sp3 "<<(st4[0]-st40[0])<<", "<<(rot[0]-rot0[0])<<", "<<(st3[0]-st30[0])<<", "<<(st3pyr[0]-st3pyr0[0])<<", "<<(st4pyr[0]-st4pyr0[0]);
fi<<" sn1 "<<(st4[1]-st40[1])<<", "<<(rot[1]-rot0[1])<<", "<<(st3[1]-st30[1])<<", "<<(st3pyr[1]-st3pyr0[1])<<", "<<(st4pyr[1]-st4pyr0[1]);
fi<<" fmn "<<(st4[2]-st40[2])<<", "<<(rot[2]-rot0[2])<<", "<<(st3[2]-st30[2])<<", "<<(st3pyr[2]-st3pyr0[2])<<", "<<(st4pyr[2]-st4pyr0[2]);
fi<<" n2 "<<(st4[3]-st40[3])<<", "<<(rot[3]-rot0[3])<<", "<<(st3[3]-st30[3])<<", "<<(st3pyr[3]-st3pyr0[3])<<", "<<(st4pyr[3]-st4pyr0[3]);
fi<<endl;
	horse.nv.setval(ip[i], a);	i++;
 } catch( char const* str ){cout << "exception: "<< str <<endl;	horse.nv.st2nv(nv1);}
	}
	
}*/
