//---------------------------------------------------------------------------
#ifndef labH
#define labH
#include <fstream>
#include <cmath>
#include <bitset>  
#include <vector>  
namespace Label {

class Metab{
protected:
      const int len;
public:
	double *isot, *disot;
	 double cont;
	int ny;
  int getlen() const {return len;}
  double total(int ini=0){
   cont=0.;
   for(int i=ini;i<len;i++) cont += isot[i];
   return cont;}

int getny() const {return ny;}
void setny(int iny) {ny=iny;}
	Metab(int l):len(l){} //
	virtual ~Metab() {}
};

class cI:public Metab{
  int nfmn2,nfmn5,n62,nq;
	public:
  double e6;
  void setnfmn2(int infs,int inq) {nfmn2=infs; nfmn5=8; n62=2; nq=inq;}
  int getnfmn2() { return nfmn2;}
void read(std::ifstream& fi) { int j(0);
  for(int i=0;i<len;i++) fi >> isot[i];  
/*      if(isot[i]>1e-7) {std::cout<<j<<") "<<i<<": "<<isot[i]<<'\n'; j++;}} */
	}
void write(std::ofstream& fi) const {
   for(int i=0;i<len;i++) fi << isot[i] << "\n";
   }
  
double fmnred(const double kf,const double kr,const double nadh,const double nad,double pfmn[],double pfmnh[]) {
 double x, sum=0; int ox,red;
 for(int k=0;k<nq;k++)	
  for(int j=0;j<n62;j++)
   for(int i=0;i<(nfmn5-2);i++) { ox=k*nfmn2+j*nfmn5+i; red=ox+2;
    if(len-red) { x=kf*pfmn[i]*nadh*isot[ox]-kr*pfmnh[i+2]*isot[red]*nad;
		disot[red] +=x; disot[ox] -=x; sum +=x; }
    else {x=kf*pfmn[i]*nadh*isot[ox]-kr*pfmnh[7]*e6*nad; disot[ox] -=x; sum +=x; }
             }
                  return sum; }
double getfs(double pfs[]) {
  double sum=0.; int j=0;
  for(int k=0;k<nq;k++)	
   for(int j=0;j<n62;j++)
    for(int i=0;i<(nfmn5);i++) { int cfs=k*nfmn2+j*nfmn5+i;
     sum +=pfs[i]*isot[cfs]; }
  return sum;}

double n2q(const double kf1,const double kr1,const double kf2,const double kr2,double pn2red[]) {
		double xf, xr, x, sum=0, kf=kf1, kr=kr1; int isub,iprod;
 for(int k=0;k<(nq-1);k++) { if(k) { kf=kf2; kr=kr2; }
  for(int j=1;j<n62;j++) 		
   for(int i=0;i<nfmn5;i++) { isub=k*nfmn2+j*nfmn5+i; iprod=isub+nfmn2-nfmn5;
        xf=kf*pn2red[j]*isot[isub]; xr=kr*(1.-pn2red[j-1])*isot[iprod]; x=xf-xr;
		disot[iprod] +=x; disot[isub] -=x; sum +=x;
		}}
		return sum; }

double n562(const double kf1,const double kr1,double pn5red[],double pn6ar[]) {
		double x, sum=0; int isub,iprod;
 for(int k=0;k<nq;k++)
  for(int j=0;j<(n62-1);j++) 		
   for(int i=1;i<nfmn5;i++) { isub=k*nfmn2+j*nfmn5+i; iprod=isub+7;
        x=kf1*pn5red[i]*(1.-pn6ar[j])*isot[isub]-kr1*(1.-pn5red[i-1])*pn6ar[j+1]*isot[iprod];
		disot[iprod] +=x; disot[isub] -=x; sum +=x;
		}
		return sum; }

double getsq(){ double sum=0.;
  for(int i=nfmn2;i<2*nfmn2;i++) sum += isot[nfmn2];
  return sum;}

double qhdiss1(cI& cr,double& qh,const double kf,const double kr) {
		double x, sum=0.; int i1;
   for(int i=0;i<(nfmn2-1);i++){ i1=2*nfmn2+i;
     x = kf*isot[i1] - kr*qh*cr.isot[i];
	 disot[i1] -=x;   cr.disot[i] += x; sum += x; }
     x = kf*isot[3*nfmn2-1] - kr*qh*cr.e6; disot[3*nfmn2-1] -=x; sum += x; 
	return sum;
	}
double qbind1(cI& qhc1,double q, const double kf,const double kr) {
		double x,sum=0.;
   for(int i=0;i<nfmn2;i++){ 
    x = kf*isot[i]*q - kr*qhc1.isot[i];
      disot[i] -=x;   qhc1.disot[i] += x; sum += x; }
   x = kf*e6*q - kr*qhc1.isot[nfmn2-1];  qhc1.disot[nfmn2-1] += x; sum += x;
	return sum;
	}	
		cI(int l):Metab(l){}
		~cI(){}
};

class cII:public Metab{
  int nfsb, nfad, nq;
	public:
  void setnfsb(int inq) {nfsb=4; nfad=3; nq=inq;}

  void read(std::ifstream& fi) {
    for(int i=0;i<len;i++) fi >> isot[i];
    //{ std::cout<<isot[i]<<' ';}std::cout<<'\n';
     }

  void write(std::ofstream& fi) const {
    for(int i=0;i<len;i++) fi << isot[i] << "\n"; }

  double fadred(const double kf,const double kr,const double suc,const double fum,double cr11) {
   double x, sum=0; int ox,red;
   for(int j=0;j<nq;j++)	
     for(int i=0;i<(nfsb);i++) { ox=j*nfsb*nfad+i*nfad; red=ox+2;
       if(red<len){x=kf*suc*isot[ox]-kr*isot[red]*fum; disot[red] +=x;}
       else {x=kf*suc*isot[ox]-kr*cr11*fum;}
       disot[ox] -=x; sum +=x;     }
   return sum; }

  double fadfs(const double kf,const double kr, double f2sr[]) {
   double x, sum(0); int ox,red, sub, pr;
   for(int j=0;j<nq;j++)	
     for(int i=0;i<(nfsb-1);i++)
       for(int k=1;k<nfad;k++) { sub=j*nfsb*nfad+i*nfad+k; pr=sub+nfad-1;
       x=kf*(1-f2sr[i])*isot[sub]-kr*isot[pr]*f2sr[i];
       disot[pr] +=x; disot[sub] -=x; sum +=x;     }
   return sum; }

  double fsq(const double kf1,const double kr1,const double kf2,const double kr2,double pbred[]) {
   double xf, xr, x, sum(0), kf(kf1), kr=(kr1); int isub,iprod;
   for(int j=0;j<(nq-1);j++) { if(j) { kf=kf2; kr=kr2; }
    for(int i=1;i<nfsb;i++) 
     for(int k=0;k<nfad;k++){
        isub=j*nfsb*nfad+i*nfad+k; iprod=isub+nfsb*nfad-nfad;
        xf=kf*pbred[i]*isot[isub];
        xr=kr*(1.-pbred[i-1])*isot[iprod]; x=xf-xr;
	disot[iprod] +=x; disot[isub] -=x; sum +=x;
		}}
   return sum; }

  double getfs() {
    double sum=0.; int j=0;
      for(int j=0;j<nq;j++)	
        for(int i=0;i<(nfsb);i++) { int cfs=j*nfsb*nfad+i*nfad+1;
          sum +=isot[cfs]; }
    return sum;}

  double fadros(const double kf) {
   double x, sum(0); int sub, pr;
   for(int j=0;j<nq;j++)	
     for(int i=0;i<(nfsb);i++)
       for(int k=1;k<nfad;k++) { sub=j*nfsb*nfad+i*nfad+1; pr=sub-1;
       x=kf*isot[sub];
       disot[pr] +=x; disot[sub] -=x; sum +=x;     }
   return sum; }

  double getsq(){
    double sum=0.;
    for(int i=nfsb*nfad;i<2*nfsb*nfad;i++) sum += isot[i];
    return sum;}

  double sqros(const double kf){
    double x, sum=0.; int sub, pr;
    for(int i=nfsb*nfad;i<2*nfsb*nfad;i++){
       sub=i; pr=sub-nfsb*nfad;
       x=kf*isot[sub];
       disot[pr] +=x; disot[sub] -=x; sum +=x;} 
    return sum;}

double qhdiss1(cII& cr,double& qh,const double kf,const double kr,double cr11) {
   double x, sum=0.; int i1;
   for(int i=0;i<cr.len;i++){ i1=2*nfad*nfsb+i;
     x = kf*isot[i1] - kr*qh*cr.isot[i];
     disot[i1] -=x;   cr.disot[i] += x; sum += x; }
     x = kf*isot[35] - kr*qh*cr11;
     disot[35] -=x;   sum += x;
   return sum;	}
double qbind1(cII& qhc1,double q, const double kf,const double kr,double cr11) {
   double x,sum=0.;
   for(int i=0;i<len;i++){
    x = kf*isot[i]*q - kr*qhc1.isot[i];
    disot[i] -=x;   qhc1.disot[i] += x; sum += x; }
    x = kf*cr11*q - kr*qhc1.isot[11];
    qhc1.disot[11] += x; sum += x;
   return sum;
	}	
		cII(int l):Metab(l){}
		~cII(){}
};

class cIII:public Metab{
	public:
   int ml, *map;
double total(int ini=0){
   cont=0.;
   for(int i=ini;i<ml;i++) cont += isot[i];
   return cont;}
int fndmap(int imap) { int i;
  for(i=0;i<ml;i++) 
    if(!(map[i]-imap)) break;
     return i;
     }
void write(std::ofstream& fi) const {
	for(int ii=0;ii<ml;ii++) fi << isot[ii] << "\n";
	}
void read(std::ifstream& fi) { int i;
    cont=0.;
    for(i=0;i<ml;i++) {fi >> isot[i]; cont += isot[i];}
	}
double percent(const int form1,const int form2) {
   double calc=0.; 
    for(int nc=0;nc<ml;nc++) if((map[nc]&form1)==form2) calc += isot[nc];
  return calc;
	}
int chkrl(int ns) {int neq=0;
   int i2=(2<<ns),i3=(3<<ns); 
    map=new int[len];
	for(int i=0;i<len;i++) 
	  if(((i&i3)-i2)&&((i&3)-2)) {map[neq]=i; neq++;} 
	 return neq;}
int chk0() {int neq=0;
    map=new int[len];
	for(int i=0;i<len;i++) 
	  if((i&3)-2) {map[neq]=i; neq++;} 
	 return neq;}	
int chkl(int ns) {int neq=0;
   int i2=(2<<ns),i3=(3<<ns); 
    map=new int[len];
	for(int i=0;i<len;i++) 
	  if((i&i3)-i2) {map[neq]=i; neq++;} 
	 return neq;}	

double c1c(int ns,const double kf,const double o2,bool fl=0) {
   double x,sum=0.; int i2;
     ns--;  int k = (1<<ns);
   for(int i=fndmap(k);i<ml;i++) 
     if (!((map[i]&k)-k)) { i2 = fndmap(map[i]-k);
	x = kf*isot[i]*o2/(0.005+o2);
	disot[i] -= x;  disot[i2] += x; 
	sum += x;}
   if(fl) {disot[ml-2] += x; sum += x;}
  return sum;}

double shift1(int ns,const double kf,const double kr) {
   double x,sum(0.); int i2;
    ns--; int i1 = (1<<ns); int i3 = (3<<ns);
     for(int i=fndmap(i1);i<ml;i++) 
        if((map[i]&i3)==i1) {i2=fndmap(map[i]+i1);
	x = kf*isot[i] -kr*isot[i2];
	disot[i] -=x;   disot[i2] += x; sum += x;}
	return sum;}	
double bhqn(int ns,const double kf1,const double kr1,const double kf2,const double kr2) {
  double x,sum=0.; int i2;
   ns--; int i1=(1<<ns), i3=(3<<ns), i5=(5<<ns), i7=(7<<ns), ij=(1<<(ns+2));
     for(int i=fndmap(i1);i<ml;i++)
      if (((map[i]&i7)==i1)) {i2=fndmap(map[i]+i1);
	x = kf1*isot[i] -kr1*isot[i2];
          disot[i] -=x;   disot[i2] += x;
	}
	else if ((map[i]&i7)==i3) {i2=fndmap(map[i]-i1+ij);
          x = kf2*isot[i] -kr2*isot[i2];
            disot[i] -=x;   disot[i2] += x; sum += x;}
   return sum;
	}	 
double qphbind(cIII& qp,double& qh, const double kf,const double kr,double bc15=0.,bool fl=0) {
	double x, sum=0.;  
	int iqp;
	for(int i=0;i<ml;i++) {	  iqp=qp.fndmap((this->map[i]<<2)+3); 
	  x = kf*qh*isot[i] -kr*qp.isot[iqp];
	  disot[i] -=x;   qp.disot[iqp] += x; sum -= x;
	}
	if(fl) {x = kf*qh*bc15 -kr*qp.isot[qp.ml-1];
	qp.disot[qp.ml-1] += x; sum -= x;}
	return sum;
	}	
double bindqn(cIII& qn,double& q, const double kf,const double kr,double bc15=0.,bool fl=0) {
		double x, sum=0.;
	for(int i=0;i<ml;i++) {	 
	  x = kf*q*isot[i] -kr*qn.isot[i];
	  disot[i] -=x;   qn.disot[i] += x;  sum -= x;//q.con -= x; 
	}
	if(fl) {x = kf*q*bc15 -kr*qn.isot[this->ml];
	qn.disot[this->ml] += x; sum -= x;}
	return sum;
	}	
double dissqn(cIII& bc1,double& qh, const double kf,const double kr,double bc15=0.,bool fl=0) {
		double x, sum=0.; 
		int i1=((len>>2)|(len>>1)); int i2;
	for(int i=0;i<bc1.ml;i++) { i2=fndmap(bc1.map[i]+i1);
	x = kf*isot[i2] -kr*qh*bc1.isot[i];
	disot[i2] -=x;   bc1.disot[i] += x; sum += x; // qh.con += x
	}
	if(fl) {x = kf*isot[this->ml-1] -kr*qh*bc15;
	disot[this->ml-1] -= x; sum += x;}
	return sum;
	}	
double qpdiss(cIII& bc1,double& q, const double kf,const double kr) {
		int iqp; double x, sum=0.;
		for(int i=0;i<ml;i++) if(!(map[i]&3)){
		iqp= bc1.fndmap(this->map[i]>>2); 
		x = kf*(isot[i]) -kr*q*bc1.isot[iqp];
		disot[i] -=x; sum += x;  
		if(iqp<bc1.ml) bc1.disot[iqp] += x;// q.con += x;
		 }
	return sum;
	}	
double shiftFeS(const double kf,const double kr,const double kfb,const double krb,const double kros,const double o2,double& so, double xc3) { 
  int i2; double x, sum=0.;
   for(int i=0;i<ml;i++)
     if ((map[i]&7)==3) { i2=fndmap(map[i]+2);
	x = kf*isot[i] -kr*isot[i2];
	disot[i] -=x;   disot[i2] += x;
	}//shift in qh2-bound from pos.2 to pos.3(FeS)
      else if ((map[i]&19)==1) { i2=fndmap(map[i]+15);
	x = kfb*isot[i] -krb*isot[i2];
	disot[i] -=x;   disot[i2] += x;	sum += x;
	}//shift in qh.-bound from pos.1 to pos.5(cytBl)
      else if ((map[i]&19)==17) { i2=fndmap(map[i]-1);
	x = kros*o2*isot[i];
	disot[i] -=x;   disot[i2] += x;	so += x*xc3;
	}//reducing of molecular oxygen
   return sum;
}
double bypass(const double kf) {
  int i, i2, k=17; double x, sum=0.;
   for(i=fndmap(k);i<ml;i++)
    if((map[i]&19)==17) { i2=fndmap(map[i]-14);
	x = kf*isot[i];
	disot[i] -=x;   disot[i2] += x; sum += x;
	}//reducing of molecular oxygen
       return sum;
}

		cIII(int l):Metab(l){}
		~cIII(){delete[] map;}
};
class Parray{
   void wstorefl (const char fn1[],int numpar,const double** m);
   void readst();
   std::string fln[22];
   double fstore[22];
public:
   int i99,i95,i90,i68;
   static double frw[];
    static int par[];
    static double nv[];
    static std::string pname[];
    double flx[22];
        Parray(){} 
        virtual ~Parray(){}
      void read(const char *fn="nv.txt",bool flg=0);
        void write (int& ifn,double,bool flg=0) const ;
	void st2nv(double st[]);
	void nv2st(double st[]){for(int i=0;i<nNV;i++) st[i] = nv[i];}
	void stat(const int NP );
	void stflx();
        double changeval (const int np,double fact){ 
         double a=nv[np]; nv[np] *= fact;
            if(np<rct) frw[np] = nv[np] * nv[rct];
		return a;}
        double setval (const int np,double a){ 
         double b=nv[np];
		nv[np] = a;
		if(np<rct) frw[np] = a * nv[rct];
		return b;
		}
	double getval(const int np) {return nv[np];}
};


class Ldistr {
        cIII pBC1q;
        cIII BC1;
        cIII qhpBC1;
        cIII BC1qn;
        cI coreI;
        cI cIq;
        cII coreII;
        cII cIIq;
       double cr11, bc15,qq, hfi, hfo,sp1,sn1,sp3,sn3, nadh,nadhc,fsq,ros,sq1, adp, ho,buf;
       double fmnh[8],fs[8],fmn[8], n5red[8], n6ar[4], n2red[4];
       double f2s[5],br[5];
       double fesr[3],c1r[3], bhr[4],sqi[4],qhi[4];
       const double c1t,c2t, c3t,vol, fc, frt, tan, tnad, tnadc,qt;
       double tshift,tex[599],ex1[2][599],sdev[599];
       int nmet, iter,ifin1,ifin2;
       double *dconc, *conc;
       vector<std::string> xname;
       double leak(double kf){return kf*exp(frt*conc[npsi])*(ho - conc[nhi]);}
       double MM(double vm,double km,double s){return (vm*s/(km+s));}
public:
       double pal,pa,pv,vo2b;
       Parray nv;
       std::string kkin, kont;
       void setdisot(double *pyinit);
       void setisot(double *pyinit);
       void chast(double *py,double tint);
       double getc3ros(){return pBC1q.percent(3,1);}
       double getc2ros(){return (coreII.getfs() + cIIq.getfs());}
int ddisolve(int istart,const double tfin, double *py,std::ostringstream& fkin);
double fiout(double t, std::ostringstream& fi,int ii=1);
double o2deliv(double pm);
double rotenone(double *py,double tint);
double antimycin(double *py);
float efit(double *py);
void ions(int nci,double P);
double jglu0();
void glufl();
void atpsyn(double katp);
void atpase(double katp);
void NaKatpase(double katp);
double gluout(double kio,double cglu){ return kio*(cglu-conc[ngluo]);}
void glycolysis();
void shutl();
	void tout(int, int&, std::ofstream&);
        int setny();
       void read (std::string fn="init");
       void seteq( double *py,double *pdydt);
       void c3calc( double *py,double *pdydt);
       void tca( double *py,double *pdydt);
       void c1calc( double *py,double *pdydt);
       void c2calc( double *py,double *pdydt);
       int readexp(int col,std::string fn="png/nadh");
       void write (std::string fn) const ;
       void distr( double *py,double *pdydt);
       double gettex(int i){return tex[i];}
	void ptp();
       double getsp3(){return sp3;}
       float dev(int kex) { float  otkl=0., a;
           for (int i=ifin1;i<ifin2;i++){ a=sdev[i]-ex1[kex][i];  otkl += a*a;}  return otkl;}
       double getconc(int i){return conc[i];}
       double setconc (int np,double a){ double b=conc[np];
		conc[np] = a; return b;}
       void setarrays();
       double getros(){return ros;}
    Ldistr():pBC1q(256), qhpBC1(64), BC1qn(64), BC1(15),coreI(15), cIq(48), c1t(0.204), coreII(11), cIIq(36), c2t(0.4), c3t(0.444), qt(4.1), vol(100.), fc(500.), frt(0.039), tan(3.1), tnad(17.0), tnadc(17.0){}
     ~Ldistr(void) {}
};
}
extern Label::Ldistr horse;
//---------------------------------------------------------------------------
#endif
 
