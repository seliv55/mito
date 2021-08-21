const int nqh=0, npsi=nqh+1, nnad=npsi+1, npyr=nnad+1, nsuc=npyr+1, nfum=nsuc+1, noaa=nfum+1, nakgm=noaa+1, nglu=nakgm+1, ngluo=nglu+1, naspm=ngluo+1, nAtp=naspm+1, neo=nAtp+1, nlac=neo+1;

void Ldistr::distr( x,dx,p1) {
double xk_i= 144.53, xNa_i= 0.01875, psio=86.266, fc=500., frt=0.039, hi=5.4e-05, ho=0.000114, krct= 1000, tan=20., tnad=20.0,
xlaco = 0.3,
xglu_o = p1,
xNa_o = 140,
xk_o = 5,
klk = 0,
 katsyn = 0.001*krct, 
 katpase = 1.9e-04*krct, 
 kGl = 0.001*krct, 
 kldh = 0.001*krct, 
 klacdif = 0.0001*krct, 
 kcI = 0.0001*krct,	
 ksdh = 0.004*krct, 
 kcIII = 0.004*krct, 
 kcs = 0.07*krct,
 kmalic = 0.0005*krct, 
 kakgsuc = 0.01*krct, 	
 kmdh = 0.9*krct, 	
 kasp_atf = 0.45*krct, 
 kasp_atr = 0.1*krct, 
 kaspout = 0.03*krct, 	
 kglu_akg = 0.0001*krct, 
 kglu_tr = 0.3*krct, 
 kgluout = 0.075*krct;
	qq=3.75-x[nqh];
	nadh= 1. - x[nnad];
	adp= tan - x[nAtp];
// TCA cycle from oaa to akg
 double v1 = kcs*x[noaa]*x[npyr]*x[nnad];
// Malic enzyme fum->pyr
 double v2 = kmalic*x[nnad]*x[nfum];
// akg to succinate
 double v3 = kakgsuc*x[nnad]*x[nakgm];
// glutamate to akg (Glutamate dehidrogenase):
 double v4 = kglu_akg*x[nglu]*x[nnad];
//leak
 double v5 = 0.005*x[npsi]*klk*exp(frt*x[npsi])*(ho - hi);
// complex I
 double v6 = kcI*nadh*tnad*qq* exp(-0.2*frt*x[npsi]);
// succinate dehydrogenase, complex II
 double v7 = ksdh*qq*x[nsuc]*(1/(1+x[noaa]/0.0002));
// complex III
 double v8 = kcIII*qq*x[nqh]* exp(-0.2*frt*x[npsi]);
// glutamate transport
 double a=125*x[ngluo]*xNa_o*xNa_o*xNa_o;
 double Tg=a*x[neo]/(a+(12500*x[ngluo]+572)*xNa_o* xNa_o + 5500*xNa_o + 110*xk_o + 1375000);
 double kt0=kglu_tr*(1/(1+x[nglu]/0.1));
 double kr0=kglu_tr*(1/(1+x[nglu]/0.1));
 double kt=kt0*exp(2*psio*frt/2), kr=kr0*exp(-psio*frt/2);
 double v9=kt*Tg;
 double nei=1.-x[neo];
	a=50*x[nglu]*xNa_i*xNa_i*xNa_i;
 double v10=kr*103180*xk_i*nei/(a+(2340*x[nglu]+10318) *xNa_i*xNa_i +25795*xNa_i+103180*xk_i+103180);
// ATP synthase
 double v11=katsyn*adp*x[npsi];
// MDH
 double v12 = kmdh*(x[nnad]*x[nfum]-nadh*x[noaa]);//
// aspartat aminotransferase: oaa + glum <-> aspm + akg
 double v13 = kasp_atf*x[noaa]*x[nglu] - kasp_atr*x[naspm]*x[nakgm];
// aspartat efflux asp ->  
 double v14= kaspout*x[naspm];
// glucose to pyruvate
 double v15= kGl*(1-x[npyr])*x[nnad];
// pyruvate to lactate
 double v16= kldh*(x[npyr]*nadh - x[nlac]*x[nnad]);
// lactate efflux/uptake
 double v17= klacdif*(xlaco - x[nlac]);
// atpase
 double v18=katpase*x[nAtp];
 double v19=kgluout*(xglu_o-x[ngluo]);
	dx[noaa] = v12*tnad - v1*tnad - v13;
	dx[npyr] = v15*tnad + v2 - v16*tnad - v1*tnad;
	dx[nakgm] = v1*tnad + v13 + v4 - v3*tnad;
	dx[nnad] = v6 - 2.*v1 - v2 - v3 - v4 - v12 - v15 + v16;
	dx[nfum] = v7 - v2 - v12*tnad;
	dx[nsuc] = v3*tnad - v7;  
	dx[nglu] = v9 - v4 - v13;
	dx[ngluo] = v19 - v9;
	dx[npsi] = 8.*v6*fc - 2.*v5*fc + 4.*v8*fc - 8.*v11*fc;
	dx[nqh] = v6 + v7 - v8;
	dx[neo] = v10 - v9;
	dx[nAtp] = v11 - v18;
	dx[naspm] = v13 - v14;
	dx[nlac] = v16*tnad + v17;
}   

