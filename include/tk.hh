//---------------------------------------------------------------------------

#ifndef TKH
#define TKH
#include "par.h"
namespace parametry {
class Doub3 {
public:
       double a[3];
       Doub3(){a[0]=0.;a[1]=0.;a[2]=0.;}
       ~Doub3(){}
};

class Doub5 {
public:
       double a[5];
       Doub5(void){a[0]=0.;a[1]=0.;a[2]=0.;a[3]=0.;a[4]=0.;}
       ~Doub5(){}
};
class TK{
       Doub5 kDe[8];
       double k[11],k0[11],v[11],v0[11];
       double e0;
public:
       void setk(Par* );
       void setree();
       void Graf(double* ) const;
       void st1fl(double* const fl) ;
       void st2fl(double* const fl);
       TK(){}
       ~TK(){}
};
class TA{
       Doub3 kDe[7];
       double k[9],k0[9],v[9],v0[9];
       double e0;
public:
       void setk( Par* );
       void setree();
       void Graf(double* ) const;
       void st1fl(double* const fl) ;
       void st2fl(double* const fl);
       TA(){}
       ~TA(){}
};
}
extern double xx[];
extern double* const px;
extern double  &r5i, &x5i, &g6i, &f6i, &g3i, &dhi, &glu, &frui,&glyci, &gtmi,&laci,&s7i, &e4i,
       &pepi, &pyri, &oaai,&citi,&coai,&mali,&akgi,&glxi;
//---------------------------------------------------------------------------
#endif
