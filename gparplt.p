#set size 0.65,0.65;
#set xrange [0.39 to 0.49];
set t pngcairo enhanced dashed font "arial,10" size 900,61.8*3;
set output "cont.png";
set multiplot layout 1,3
set style line 1 lt 1 lw 1.5 lc rgb "#3465a4"
set style line 2 lt 1 lw 1.5 lc rgb "#f57900"
set xlabel "Succinate (relative)" offset 0,0.7;
set xtics nomirror;
set ytics nomirror;
set border 3;
eflow=3; sqo=eflow+1; fsc2=sqo+1; qsc2=fsc2+1; fum=qsc2+1; bypas=fum+1; qh=bypas+1; psi=qh+1; fc1=psi+1; rosc2=fc1+1; rosc3=rosc2+1; fc2=rosc3+1;
set ylabel "C3 ROS rate (relative)" offset 2;
set ytics auto;
plot '00000' using ($1*10):sqo w l ls 2 not;
set ylabel "C2 ROS rate (relative)" offset 2;
plot '00000' using ($1*10):rosc3 w l ls 2 not;
set ylabel "ROS_c2" offset 2;
plot '00000' using ($1*10):rosc2 w l ls 2 not;

