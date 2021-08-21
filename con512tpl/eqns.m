
function dx=eqns( x,p1) 
 xk_i= 144.53; xNa_i= 0.01875; psio=86.266; fc=500.; frt=0.039; hi=5.4e-05; ho=0.000114; krct= 1000; tan=20.; tnad=20.0;
xlaco = 0.3;
xglu_o = p1;
xNa_o = 140;
xk_o = 5;
klk = 0;
 katsyn = 0.001*krct; 
 katpase = 1.9e-04*krct; 
 kGl = 0.001*krct; 
 kldh = 0.001*krct; 
 klacdif = 0.0001*krct; 
 kcI = 0.0001*krct;	
 ksdh = 0.004*krct; 
 kcIII = 0.004*krct; 
 kcs = 0.07*krct;
 kmalic = 0.0005*krct; 
 kakgsuc = 0.01*krct; 	
 kmdh = 0.9*krct; 	
 kasp_atf = 0.45*krct; 
 kasp_atr = 0.1*krct; 
 kaspout = 0.03*krct; 	
 kglu_akg = 0.0001*krct; 
 kglu_tr = 0.3*krct; 
 kgluout = 0.075*krct;
	qq=3.75-x(1);
	nadh= 1. - x(3);
	adp= tan - x(12);
% TCA cycle from oaa to akg
  v1 = kcs*x(7)*x(4)*x(3);
% Malic enzyme fum->pyr
  v2 = kmalic*x(3)*x(6);
% akg to succinate
  v3 = kakgsuc*x(3)*x(8);
% glutamate to akg (Glutamate dehidrogenase):
  v4 = kglu_akg*x(9)*x(3);
%leak
  v5 = 0.005*x(2)*klk*exp(frt*x(2))*(ho - hi);
% complex I
  v6 = kcI*nadh*tnad*qq* exp(-0.2*frt*x(2));
% succinate dehydrogenase, complex II
  v7 = ksdh*qq*x(5)*(1/(1+x(7)/0.0002));
% complex III
  v8 = kcIII*qq*x(1)* exp(-0.2*frt*x(2));
% glutamate transport
  a=125*x(10)*xNa_o*xNa_o*xNa_o;
  Tg=a*x(13)/(a+(12500*x(10)+572)*xNa_o* xNa_o + 5500*xNa_o + 110*xk_o + 1375000);
  kt0=kglu_tr*(1/(1+x(9)/0.1));
  kr0=kglu_tr*(1/(1+x(9)/0.1));
  kt=kt0*exp(2*psio*frt/2), kr=kr0*exp(-psio*frt/2);
  v9=kt*Tg;
  nei=1.-x(13);
	a=50*x(9)*xNa_i*xNa_i*xNa_i;
  v10=kr*103180*xk_i*nei/(a+(2340*x(9)+10318) *xNa_i*xNa_i +25795*xNa_i+103180*xk_i+103180);
% ATP synthase
  v11=katsyn*adp*x(2);
% MDH
  v12 = kmdh*(x(3)*x(6)-nadh*x(7));%
% aspartat aminotransferase: oaa + glum <-> aspm + akg
  v13 = kasp_atf*x(7)*x(9) - kasp_atr*x(11)*x(8);
% aspartat efflux asp ->  
  v14= kaspout*x(11);
% glucose to pyruvate
  v15= kGl*(1-x(4))*x(3);
% pyruvate to lactate
  v16= kldh*(x(4)*nadh - x(14)*x(3));
% lactate efflux/uptake
  v17= klacdif*(xlaco - x(14));
% atpase
  v18=katpase*x(12);
  v19=kgluout*(xglu_o-x(10));
	dx(1) = v6 + v7 - v8;
	dx(2) = 8.*v6*fc - 2.*v5*fc + 4.*v8*fc - 8.*v11*fc;
	dx(3) = v6 - 2.*v1 - v2 - v3 - v4 - v12 - v15 + v16;
	dx(4) = v15*tnad + v2 - v16*tnad - v1*tnad;
	dx(5) = v3*tnad - v7;  
	dx(6) = v7 - v2 - v12*tnad;
	dx(7) = v12*tnad - v1*tnad - v13;
	dx(8) = v1*tnad + v13 + v4 - v3*tnad;
	dx(9) = v9 - v4 - v13;
	dx(10) = v19 - v9;
	dx(11) = v13 - v14;
	dx(12) = v11 - v18;
	dx(13) = v10 - v9;
	dx(14) = v16*tnad + v17;
%endfunction
end
