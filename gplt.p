set t pngcairo enhanced dashed font "arial,10" size 1000,400;
set xlabel "Time (min)" offset 0,0.7;
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
#//set ylabel "nadh (nmol/mg)" offset 2;
#//plot '011110rbm' u ($1):($2) every 3::3::293 ps 0.25 lc rgb "#3465a4" not,'011110rbm' u ($1):($3) every 3::3::293 ps 0.25 lc rgb "#f57900" not;
set ytics autofreq
set ylabel "ROS rate in CII (relative)" offset 2;
plot fn1 using ($1):fsc2 w l ls 1 not, fn2 using ($1):fsc2 w l ls 2 not;
set ylabel "ΔΨ (mV)" offset 2;
plot fn1 using ($1):psi w l ls 1 not, fn2 using ($1):psi w l ls 2 not;
set ylabel "e-flow (1/s)" offset 2;
plot fn1 using ($1):eflow w l ls 1 not, fn2 using ($1):eflow w l ls 2 not;
set ylabel "QH_2 (nmol/mg)" offset 2;
plot fn1 using ($1):qh w l ls 1 not, fn2 using ($1):qh w l ls 2 not;
set ylabel "ROS rate in CIII (relative)" offset 2;
plot fn1 using ($1):sqo w l ls 1 not, fn2 using ($1):sqo w l ls 2 not;
#set ylabel "[Ca²⁺]μM" offset 2;
set ylabel "ROS in C2" offset 2;
plot fn1 using ($1):rosc2 w l ls 1 not, fn2 using ($1):rosc2 w l ls 2 not;
set ylabel "ROS in C3" offset 2;
#set yrange [-0.1:0.1]
plot fn1 using ($1):rosc3 w l ls 1 not, fn2 using ($1):rosc3 w l ls 2 not;
set ylabel "qsc2" offset 2;
set yrange [*:*]
plot fn1 using ($1):qsc2 w l ls 1 not, fn2 using ($1):qsc2 w l ls 2 not;
unset multiplot
set output "Ca.png";
set t pngcairo enhanced dashed font "arial,10" size 400,400;
set multiplot layout 2,1;
set yrange [0:220]
set ylabel "ΔΨ (mV)" offset 2;
plot fn1 using ($1):psi w l ls 1 not, fn2 using ($1):psi w l ls 2 not;
set yrange [0:20]
set ylabel "[Ca²⁺]μM" offset 2;
plot fn1 using ($1):ca w l ls 1 not, fn2 using ($1):ca w l ls 2 not;


