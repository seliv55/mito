#set size 0.65,0.65;
set t pngcairo enhanced dashed font "arial,12";
set output "GluTr.png";
set style line 1 lt 1 lw 1.5 lc rgb "#3465a4"
set style line 2 lt 1 lw 1.5 lc rgb "#f57900"
set xlabel "Glutamate, mM" offset 0,0.7;
#set xrange [0:200]
set xtics nomirror;
set ytics nomirror;
set border 3;
#eflow=3; sqo=eflow+1; fsc2=sqo+1; qsc2=fsc2+1; fum=qsc2+1; bypas=fum+1; qh=bypas+1; psi=qh+1; fc1=psi+1; rosc2=fc1+1; rosc3=rosc2+1; fc2=rosc3+1;
set ytics auto;
set ylabel "Glu transport rate, relative" offset 2;
plot 'glutr' using ($1):2 w l ls 2 not;

