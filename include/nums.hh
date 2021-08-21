extern const int  nqh, npsi, nnad, npyr, nsuc, nfum, noaa, nakgm, nglu, ngluo, naspm, nca, nc1ros, nc2ros, nc3ros, nk_i, npsio, nNa_i, nAdp, nAtp, neo, nei, nlac, numx;
extern int kmax,kount,ifn,istor,nrhs;
extern const int  NN, qp_FS, rqp_FS, FS_c1, rFS_c1, qp_bl, rqp_bl, bl_bh, rbl_bh, bh_qn1, rbh_qn1, bh_qn2, rbh_qn2, qHbnd, rqHbnd, qnbnd, rqnbnd, qpdis, rqpdis, qhnds,rqhnds, vc1c, bypas, vros, vfadf, vfadr, vbq1, vrbq1, vbq2, vrbq2, vqdis, vrqdis, vqbind, vrqbind, vfred, vrfred, vn56, vrn56, vn2qn1, vrn2qn1, vn2qn2, vrn2qn2, vndis, vrndis, vpbind, vrpbind, ptp, vlk, vatsyn, vmalic, vmdh, vK, akgsuc, vglu_akg, vatpase, vglu_tr, vlacdif, vgluout, vcI, vsdh, vcIII, vcs, vasp_atf, vasp_atr, vGl, vldh, vaspout, vperos, rct, ntmax, katase, glu_o, laco, k_o, Na_o, nNV;
void arch(double *y,int ipar,double pint,double step);
void jsvd(double *y,Vec_DP& w,Mat_DP& u,Mat_DP& v,Mat_DP& a, double cond,const int n);
double neuton(double *y,Vec_DP& c,Vec_DP& w,Mat_DP& u,Mat_DP& v,const int n);
void correct(double *y,Vec_DP& c,Vec_DP& w,Mat_DP& u,Mat_DP& v,const int n);
void init_ydot(double *y,int ipar,double dp,Vec_DP&,double*,const int n);
void ydot(double *y,int ipar,double dp,Vec_DP&,const int n);
void jlu(double *y,Mat_DP& u);
extern const int fsdh, fsm, fc1c, fme, fpyrt, fcs, ftca, fuoa, fatps, flk, fbp, fc1, fc2, flkc1, nflx;
extern double tfac, xi1, xi2, x01, x02, x03, x04;

