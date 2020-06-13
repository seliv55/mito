//---------------------------------------------------------------------------
#ifndef labH
#define labH
#include <fstream>
#include <cmath>
namespace Label {

class Metab{
protected:
      const int len;
public:
	double *isot, *disot;
	 double cont;
	int ny;

  int getlen() const {return len;}
	Metab(int l):len(l){} //
	virtual ~Metab() {}
};

class cI:public Metab{
	public:
  void total(){
   cont=0.;
   for(int i=0;i<len;i++) cont += isot[i];
   }
void read(std::ifstream& fi) {
  for(int i=0;i<len;i++) fi >> isot[i]; 
	}
void write(std::ofstream& fi) const {
   for(int i=0;i<len;i++) fi << isot[i] << "\n";
	}
double fmnred(const double hi,const int acc,const int red,const double kf,const double kr, const double nadh,const double nad) {
		double x;
		x=kf*hi*nadh*isot[acc]-kr*isot[red]*nad;
		disot[red] +=x; disot[acc] -=x; 
			return x; 
	}
	
double redox1(int don, int acc,double kf,double kr) {
		double x=kf*isot[don]-kr*isot[acc];
		disot[acc] +=x; disot[don] -=x;
			return x;  }
double redox2(int don, double acc,double kf,double kr) {
		double x=kf*isot[don]-kr*acc;
		disot[don] -=x;
			return x;  }
			
double qpdiss1(Metab& qncom1,double& qh,int iqnp,int iqn, const double kf,const double kr) {
		double x;
	x = kf*(isot[iqnp]) -kr*qh*qncom1.isot[iqn];
	disot[iqnp] -=x;   qncom1.disot[iqn] += x; //qh.con += x; 
	return x;
	}
	
double qpbind1(Metab& qp,double& q,double qn,int iqnp, const double kf,const double kr) {
		double x;  
	x = kf*q*qn -kr*qp.isot[iqnp];
		qp.disot[iqnp] += x; 
		return x;
	}	
	
		cI(int l):Metab(l){}
		~cI(){}
};
class cIII:public Metab{
	public:
   int ml, *map;
  void total(){
   cont=0.;
   for(int i=0;i<ml;i++) cont += isot[i];
   }
int fndmap(int imap) { int i;
  for(i=0;i<ml;i++) 
    if(!(map[i]-imap)) break;
     return i;
     }
void write(std::ofstream& fi) const {
	for(int i=0;i<ml;i++) fi << isot[i] << "\n";
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
	for(int i=1;i<len;i++) 
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

double c1c(int ns,const double kf,const double o2) {
   double x,sum=0.; int i2;
     ns--;  int k = (1<<ns);
   for(int i=fndmap(k);i<ml;i++) 
     if (!((map[i]&k)-k)) { i2 = fndmap(map[i]-k);
	x = kf*isot[i]*o2/(0.005+o2);
	disot[i] -= x; if(i2) disot[i2] += x; 
	sum += x;}
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
double qphbind(cIII& qp,double& qh, const double kf,const double kr) {
		double x, sum=0.;  
		int iqp;
	for(int i=0;i<ml;i++) { iqp=qp.fndmap((map[i]<<2)+3); 
	x = kf*qh*isot[i] -kr*qp.isot[iqp];
	disot[i] -=x;   qp.disot[iqp] += x; sum -= x;
	}
	return sum;
	}	
double bindqn(cIII& qn,double& q,double& c0, const double kf,const double kr,bool flg) {
		double x, sum=0.;
		int i1;
	for(int i=1;i<ml;i++) {  i1=qn.fndmap(this->map[i]);
	x = kf*q*isot[i] -kr*qn.isot[i1];
	disot[i] -=x;   qn.disot[i1] += x;  sum -= x;//q.con -= x; 
	}
	if(flg) x = kf*q*isot[0] -kr*c0; 
	else {x = kf*q*isot[0] -kr*qn.isot[0]; qn.disot[0] += x;}
	disot[0] -= x; sum -= x;
	return sum;
	}	
double dissqn(cIII& bc1,double& qh, const double kf,const double kr) {
		double x, sum=0.; 
		int i1=(bc1.len|(len>>1)); int i2;
	for(int i=0;i<bc1.ml;i++) { i2=fndmap(bc1.map[i]+i1);
	x = kf*isot[i2] -kr*qh*bc1.isot[i];
	disot[i2] -=x;   bc1.disot[i] += x; sum += x; // qh.con += x
	}
	return sum;
	}	
double qpdiss(cIII& bc1,double& q,double& c0, const double kf,const double kr,bool flg) {
		int iqp; double x, sum=0.;
		for(int i=0;i<ml;i++) if(!(map[i]&3)){
		iqp= bc1.fndmap(this->map[i]>>2); 
		x = kf*(isot[i]) -kr*q*bc1.isot[iqp];
		disot[i] -=x;   bc1.disot[iqp] += x; sum += x;// q.con += x;
		 }
	if(flg) {x = kf*c0 -kr*q*bc1.isot[0]; bc1.disot[0] += x; sum += x;}
	return sum;
	}	
double shiftFeS(const double kf,const double kr,const double kfb,const double krb,const double kros,const double o2,double& so) { 
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
	disot[i] -=x;   disot[i2] += x;	so += x;
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
		~cIII(){}
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
    double qpqn, rqpqn,c3sq[7],c1sq[7],fmnh[7],rn2[7], flx[22];
        Parray(){} 
        virtual ~Parray(){}
      void read(const char *fn="nv.txt",bool flg=0);
        void write (int& ifn,double,bool flg=0) const ;
	void st2nv(double st[]);
	void nv2st(double st[]);
	void stat(const int NP );
	void stflx();
        double changeval (const int np,double fact){ 
         double a=nv[np];
		nv[np] *= fact;
		if(np<rct) frw[np] = nv[np] * nv[rct];
		return a;
		}
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
        cI qn1;
        cI qnp1;
       double dpdt,sp1,sn1,sp3,sn3, fmn, n2;
       const double c1t,c3t,vol, buf, fc, frt, tan, tnad,qt,ho,hi;
       double tshift,tex[599],ex1[599];
	const int in10100;
const int i1100000,i1100101,i1100011,i1110100,i1110001,i1101001,i1011100,i1011001;
       int nmet;
       double *dconc;
       double *conc;
double leak(double kf){return kf*exp(frt*conc[npsi])*(ho - hi);}
       void setdisot(double *pyinit);
       double MM(double vm,double km,double s){return (vm*s/(km+s));}
public:
       double pal,pa,pv,vo2b;
       Parray nv;
        void setisot(double *pyinit);
double fiout(double t, std::ofstream& fi,int ii, int crv=1);
double o2deliv(double pm);
	void tout(int, int&, std::ofstream&);
	void priout();
        int setny();
       void read (std::string fn="init");
       int readexp(int col,std::string fn="png/nadh");
       void write (char fn[]) const ;
       void distr( double *py,double *pdydt,double,double);
       double getdpdt()const{return dpdt;}
       double gettex(int i){return tex[i];}
       double getsp3(){return sp3;}
       double getconc(int i){return conc[i];}
       double setconc (int np,double a){ double b=conc[np];
		conc[np] = a; return b;}
 void getros(double xi[]){ xi[0] += sp3; xi[1] += sn1;  xi[2] += fmn; xi[3] += n2;}
	void wros(std::ofstream& fi){ fi<<sp3<<" "<<sn1<<" "<<fmn<<" "<<n2<<std::endl;}
	Ldistr():pBC1q(256), qhpBC1(64), BC1qn(64), BC1(16),qn1(1), qnp1(8), c1t(0.204),c3t(0.4),qt(7.1), hi(5.4e-05), ho(0.000114), vol(100.),buf(1000000.),fc(500.),frt(0.039),tan(1.),tnad(17.0),in10100(0),i1100000(0),i1100101(1),i1100011(2),i1110100(3),i1110001(4),i1101001(5),i1011100(6),i1011001(7){}
	~Ldistr(void) {}
};
}
extern Label::Ldistr horse;
//---------------------------------------------------------------------------
#endif
 
