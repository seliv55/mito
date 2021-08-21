set t pngcairo enhanced dashed font "arial,10" size 1000,400;
set xlabel "Time (s)" offset 0,0.7;
set xrange [0 to ax]
set output fno;
set multiplot layout 2,4;
set style line 1 lt 1 lw 1.5 lc rgb "#3465a4"
set style line 2 lt 1 lw 1.5 lc rgb "#f57900"
set xtics 1.0 nomirror;
set ytics 3 nomirror;
set border 3;
set key horizontal at graph 0.5,1.1 left Right samplen 1;
eflow=2; sqo=eflow+1; fsc2=sqo+1; qsc2=fsc2+1; fum=qsc2+1; bypas=fum+1; qh=bypas+1; psi=qh+1; fc1=psi+1; rosc2=fc1+1; rosc3=rosc2+1; fc2=rosc3+1; ca=fc2+1;
Glu_m=ca+1; OAAm=Glu_m+1; psio=OAAm+1; Na_i=psio+1; OAAc=Na_i+1; Glu_o=OAAc+1; Glu_i=Glu_o+1; ATP=Glu_i+1;
set ytics autofreq
set ylabel "Glu mito (mM)" offset 2;
plot fn1 using ($1):Glu_m w l ls 1 not;
set ylabel "OAA mito (mM)" offset 2;
plot fn1 using ($1):OAAm w l ls 1 not;
set ylabel "ΔΨ cell (mV)" offset 2;
plot fn1 using ($1):psio w l ls 1 not;
set ylabel "Na+ in (mM)" offset 2;
plot fn1 using ($1):Na_i w l ls 1 not;
set ylabel "OAA cyt (mM)" offset 2;
plot fn1 using ($1):OAAc w l ls 1 not;
set yrange [* to *]
set ylabel "Glu out, mM" offset 2;
plot fn1 using ($1):Glu_o w l ls 1 not;
set ylabel "Glu in, mM" offset 2;
plot fn1 using ($1):Glu_i w l ls 1 not;
set ylabel "ATP, mM" offset 2;
plot fn1 using ($1):ATP w l ls 1 not;


