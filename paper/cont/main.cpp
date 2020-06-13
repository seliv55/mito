//---------------------------------------------------------------------------
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "nr.h"
#include "StiffIntegratorT.h"
#include "NonStiffIntegratorT.h"
#include "nums.hh"
#include "modlab.h"
#include "analis.h"
#include "main.hh"
//---------------------------------------------------------------------------
//#pragma package(smart_init)


using namespace std;
using namespace Label;
extern vector<double> vx;
Ldistr horse;
//int  method=0;
const int NN=horse.setny()+numx, numex=horse.readexp(2,"a011110rbm.txt");
double ystart[500];
DP dxsav;   // defining declarations
DP *yy0, *yy1;
int kmax,kount;
Vec_DP *xp_p;
Mat_DP *yp_p;
Vec_INT *ija_p;
Vec_DP *sa_p;
int nrhs,ifn=0,istor=20001;   // counts function evaluations
Vec_DP *x_p;
Mat_DP *d_p;
    time_t ts,tf,tcal,tfirst;
double xi1, xi2,x01,x02,x03,x04,st4[4],rot[4],st3[4],st3pyr[4],st4pyr[4],pyrot[4];
    const double tfac=2100.;

 double solve(double *py,double tint){
    double xi, y[NN];
    for(int i=0;i<NN;i++) y[i]=py[i];
       horse.readexp(2,"a011110rbm.txt");//read experiment
    ts=clock();
     try { 
// xi1=tisolve(horse.nv.nv[ntmax], y,99);
 horse.chast( py,tint);
// xi=d02nefsv(horse.nv.nv[ntmax], y,numex,crv,ros);
    tf=clock()-ts;
}catch( char const* str ){cout <<"exception: "<<str<<endl; tf=1.e11; xi=1.e11;}
//    horse.priout();
//    horse.nv.stflx();
    return horse.getc3ros();}//py[265];}//horse.getconc(npsi)
    
double cont(double *py,double tint){
 const int nss(20);
 const int param[]={0,-1,1,2,3,4,5,6,7,8,9,10,11, 12, 13, 14, 15, 16, 17, 
 18, 19, 20, 37, -1};
   double a,dif,fact(1.0),f1(1.02),mx(0.),ss[2][nss],df[33], pdf[33];
   const float pmin(0.001),pint(0.3), dp=pint/(float)nss;
   stringstream fn; 
   ofstream fo;
   int chpar=suc_o;//vfadf;//pyr_o;//
   double vpar=horse.nv.nv[chpar];
   int j=0;
   double sum(0.),ss1(0.),ss2(0.);
    for(int k=0;k<2;k++){
   cout<<"nss="<<nss<<endl;
      for(int ii=0;ii<(nss+1);ii++){
        if(k) horse.nv.setval(chpar,(pmin+(pint-ii*dp)));
        else horse.nv.setval(chpar,(pmin+ii*dp));
        horse.setisot(py); 
//         horse.read("i1");// read init values
//     ss[k][ii]= horse.rotenone(py,tint);
     ss[k][ii]= solve(py,tint);
     float ft = float(tf)/CLOCKS_PER_SEC;
     cout<<ii<<" "<<horse.nv.getval(chpar)<<" ψ="<<ss[k][ii]<<" t="<<ft<<'\n';
     fn<<horse.nv.getval(chpar)<<" "<<horse.kont;//ss[k][ii]<<'\n';//
     if(ii==nss) break;
     if(ft>3){
       fo.open("nor"); fo<<horse.kkin<<endl; fo.close(); 
//       horse.rotenone(py);
//       fo.open("rot"); fo<<horse.kkin<<endl; fo.close();
ostringstream pltt;
pltt<< "gnuplot -e \"fn1='nor'; fn2='nor'; fno='kin/suc.png'; ax="<<tint<<"\" gplt.p";
       int sys=system(pltt.str().c_str());
     }
     } }
     string aaa=fn.str(); cout<<aaa;
   for(int i=0;i<nss;i++) { dif=abs(ss[0][i]-ss[1][nss-i]);
     if(dif>0.1) cout<<i<<" "<<dif<<"; ";
     sum += dif;
   } cout<< endl;
   cout<<" sum="<<sum<<endl;
     fo.open("00000"); fo<<aaa<<endl; fo.close();
return 1.;}

double comb(double *py,int wr){
 const int ndf=33, nss=30;//1;//
 const int param[]={20, -1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 37, -1};
   double fact(1.1),mx(0.),a,dif,df[ndf],pdf[ndf];
   ofstream fo("parameters");
   int chpar=suc_o;//pyr_o;
   double vpar=horse.nv.nv[chpar];
   for(;;){
   int j=0;
   while(param[j]>=0){ a = horse.nv.changeval(param[j],fact);
   double sum(0.),ss1(0.),ss2(0.); 
    for(int ii=0;ii<50;ii++){//while(ss2<0.31) {// 
      horse.nv.changeval(chpar,1.1);
       horse.setisot(py); 
         horse.read("i1");// read init values
     ss1= solve(py,5.); tcal = tf;//control
     if(ss1<-0.5) fo<<"**No steady state!**"<<endl;
/**/       horse.setisot(py); 
         horse.read("i0");// read init values
     ss2= solve(py,5.); tcal = tf;//control
     if(ss2<-0.5) fo<<"**No steady state!**"<<endl;
     dif=ss1-ss2; dif *= dif; sum += dif; 
   cout<<"p="<<horse.nv.nv[chpar]<<"; s1="<<ss1<<"; s2="<<ss2<<endl;
   fo<<"p="<<horse.nv.nv[chpar]<<"; s1="<<ss1<<"; s2="<<ss2<<endl;
   }
   df[j]=sum; pdf[j]=horse.nv.nv[param[j]];
   cout<<param[j]<<": "<<horse.nv.nv[param[j]]<<", dif="<<sum<<endl;
   fo<<param[j]<<": "<<horse.nv.nv[param[j]]<<", dif="<<sum<<endl;
   horse.nv.setval(param[j],a);    horse.nv.setval(chpar,vpar); 
      j++;
   } int imx(-1);
   for(int i=0;i<j;i++) if(mx<df[i]) {mx=df[i]; imx=i;} 
   if(imx>=0){horse.nv.setval(param[imx],pdf[imx]); 
      cout<<"change! "<<param[imx]<<": "<<pdf[imx]<<endl;
         fo<<"change! "<<param[imx]<<": "<<pdf[imx]<<endl;
           horse.nv.write(ifn,xi1);
         }
   fact=1/fact;}
return 1.;}

double Ldistr::rotenone(double *py,double tint){
 double a1=nv.setval(vn2qn1, 0.0), a2=nv.setval(vrn2qn1, 0.0), a3=nv.setval(vn2qn2, 0.0), a4=nv.setval(vrn2qn2, 0.0);
  double x=solve(py,tint);
  nv.setval(vn2qn1,a1); nv.setval(vrn2qn1,a2), nv.setval(vn2qn2, a3), nv.setval(vrn2qn2, a4); 
  return x;}

double Ldistr::antimycin(double *py){
  double a1=nv.setval(qnbnd, 0.0), a2=nv.setval(rqnbnd, 0.0),
     a3=nv.setval(qhnds, 0.0), a4=nv.setval(rqhnds, 0.0);
  double x=solve(py,5);
  nv.setval(qnbnd, a1); nv.setval(rqnbnd, a2);
     nv.setval(qhnds, a3); nv.setval(rqhnds, a4);
  return x;}
    
float Ldistr::efit(double *py){
  float xii(0);
  // solve(py); int tf1=tf;  xii=dev(0);
    rotenone(py,5.);  xii += dev(1); // tf += tf1;
     return xii;}

int main( int argc, char *argv[] ){
   cout<<"NN="<<NN<<'\n';
   if(argc==1) cout<<"run: ./a.out [parameter:'f'-fit,'s'-statistic,'a'-sensitivity]"<<endl; 
   char fn[22];
 //define last parameter-file:
   for(int i=ifn+1;;i++) { sprintf(fn,"%i",i);
	   ifstream checkfi(fn);
	   if(!checkfi.good()) { ifn=i; break;}
	   checkfi.close();
	   }
 //statistic:
   horse.setarrays();
   if ((argc==2)&&(argv[1][0]=='s'))  { horse.nv.stat(ifn); return 0; }
 //calculation: 
 //    deliv_();
   double *py = &ystart[0], tint=7;  cout.precision(3);
   horse.setisot(&ystart[0]); 
   horse.read("i1");//"i0" read init values
   horse.nv.read("1",1);//0-pyr+mal; 1-suc
   xi1=xi2=1.; 
//        double xi=comb(py);
   double x=horse.rotenone(py,tint);// solve(&ystart[0],7.);//
   ofstream fo("nor"); fo<<horse.kkin<<endl; fo.close(); 
//   horse.rotenone(py); cout<<"time="<<tf/double(CLOCKS_PER_SEC)<<endl;
   fo.open("rot"); fo<<horse.kkin<<endl; fo.close(); 
ostringstream pltt;
pltt<< "gnuplot -e \"fn1='nor'; fn2='nor'; fno='kin/suc.png'; ax="<<tint<<"\" gplt.p";
       int sys=system(pltt.str().c_str());
   horse.nv.write(istor,xi1,0);
   double xi=cont(&ystart[0],7);    //sys=system("gnuplot gparplt.p");
 //start fitting
   Analis analis; srand(time(NULL));

//
   if((argc==2)&&(argv[1][0]=='f')) analis.coord(ystart);
//   if((argc==2)&&(argv[1][0]=='a')) analis.sensitivity(ystart);
}

