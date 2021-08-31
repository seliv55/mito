#set t pngcairo enhanced dashed font "arial,10" size 900,68*6; set output 'kin/suc.png';
set t svg enhanced dashed font "arial,12" size 708,68*5; set output 'kin/suc.svg';
set xlabel "Time (s)" offset 0,0.7;
set xrange [0 to *]
set multiplot layout 2,3;
set style line 1 lt 1 lw 2 lc rgb "#204a87"
set style line 2 lt 1 lw 2 lc rgb "#4e9a06"
set style line 3 lt 1 lw 2 lc rgb "#f57900"
set style line 4 lt 1 lw 2 lc rgb "#2e3436"
set xtics 10 nomirror;
set ytics 3 nomirror;
set border 3;
set key ver at graph 0.85,0.8 left Right samplen 1;
eflow=2; sqo=eflow+1; fsc2=sqo+1; nad=fsc2+1; fum=nad+1; suc=fum+1; qh=suc+1; psi=qh+1; Vrosc1=psi+1; rosc2=Vrosc1+1; rosc3=rosc2+1; fc2=rosc3+1; ca=fc2+1;
Glu_m=ca+1; OAAm=Glu_m+1; psio=OAAm+1; Na_i=psio+1; OAAc=Na_i+1; Glu_o=OAAc+1; Glu_i=Glu_o+1; ATP=Glu_i+1; pH=ATP+1;
#//plot '011110rbm' u ($1):($2) every 3::3::293 ps 0.25 lc rgb "#3465a4" not;
#--------------------------------------------------------------
set ytics 0.2
set yrange [-0.1 to *]
set label 1 "A" at graph 0.05,1.0 font "arial,12"
set ylabel "QH_2 fraction" offset 2;
plot "ph6" using ($1):qh w l ls 3 t"1", fn2 using ($1):qh w l ls 1 t"2", fn1 using ($1):qh w l ls 2 t"3"#, "ph9" using ($1):qh w l ls 4 t"pH=9"
#plot fn1 using ($1):qh w l ls 2 t"pH=8"

set ytics 1
set yrange [5.5 to *]
set label 1 "B" at graph 0.05,1.0 font "arial,12"
set ylabel "pH" offset 2;
plot fn1 using ($1):pH w l ls 2 not, fn2 using ($1):pH w l ls 1 not, "ph6" using ($1):pH w l ls 3 not"pH=6"#, "ph9" using ($1):pH w l ls 4 not"pH=9"

set ytics 50
set yrange [-5 to *]
set label 1 "C" at graph 0.05,1.0 font "arial,12"
set ylabel "Δψ (mV)" offset 2;
plot fn1 using ($1):psi w l ls 2 not, fn2 using ($1):psi w l ls 1 not, "ph6" using ($1):psi w l ls 3 not"pH=6"#, "ph9" using ($1):psi w l ls 4 not"pH=9"

set ytics 0.05
set yrange [0 to *]
set label 1 "D" at graph 0.05,1.0 font "arial,12"
set ylabel "CI & CII ROS rate" offset 2;
plot fn1 using ($1):Vrosc1 w l ls 2 not, fn2 using ($1):Vrosc1 w l ls 1 not, "ph6" using ($1):Vrosc1 w l ls 3 not"pH=6", fn1 using ($1):fsc2 w l ls 2 dt 2 not, fn2 using ($1):fsc2 w l ls 1 dt 4 not, "ph6" using ($1):fsc2 w l ls 3 not"pH=6"#, "ph9" using ($1):Vrosc1 w l ls 4 not"pH=9"

set ytics 0.2
set yrange [-0.05 to *]
set label 1 "E" at graph 0.05,1.0 font "arial,12"
set ylabel "CIII ROS rate" offset 2;
plot fn1 using ($1):sqo w l ls 2 not, fn2 using ($1):sqo w l ls 1 not, "ph6" using ($1):sqo w l ls 3 not"pH=6"#, "ph9" using ($1):sqo w l ls 4 not"pH=9"

#set ytics 0.1
#set label 1 "F" at graph 0.05,1.0 font "arial,12"
#set ylabel "CII ROS rate" offset 2;
#plot fn1 using ($1):fsc2 w l ls 2 not, fn2 using ($1):fsc2 w l ls 1 not, "ph6" using ($1):fsc2 w l ls 3 not"pH=6", "ph9" using ($1):fsc2 w l ls 4 not"pH=9"

set ytics auto
set label 1 "F" at graph 0.05,1.0 font "arial,12"
set ylabel "Ca^2^+ (nmol/mg)" offset 2;
plot fn1 using ($1):ca w l ls 2 not, fn2 using ($1):ca w l ls 1 not, "ph6" using ($1):ca w l ls 3 not"pH=6"#, "ph9" using ($1):ca w l ls 4 not"pH=9"
#set label 1 "B" at graph 0.05,1.0 font "arial,12"
#set ylabel "NAD (relative)" offset 2;
#plot fn1 using ($1):rosc3 w l ls 2 not#, fn2 using ($1):rosc3 w l ls 3 not, "rot1" using ($1):rosc3 w l ls 1 not;
#set ytics autofreq
#set yrange [-0.0 to *]
#set label 1 "C" at graph 0.05,1.0 font "arial,12"
#set ylabel "Glu mito (μmol/g)" offset 2;
#plot fn1 using ($1):Glu_m w l ls 2 not#, fn2 using ($1):Glu_m w l ls 3 not, "rot1" using ($1):Glu_m w l ls 1 not;
#set ylabel "ΔΨ (mV)" offset 2;
#plot fn1 using ($1):psi w l ls 1 not, fn2 using ($1):psi w l ls 2 not;
#set ylabel "e-flow (1/s)" offset 2;
#plot fn1 using ($1):eflow w l ls 1 not, fn2 using ($1):eflow w l ls 2 not;

unset multiplot


