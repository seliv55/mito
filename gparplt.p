#set size 0.65,0.65;
set t pngcairo enhanced dashed font "arial,12" size 900,61.8*6;
#set t svg enhanced dashed font "arial,10" size 600,61.8*3;
set output "cont.png";
set multiplot layout 1,2
set style line 1 lt 1 lw 1.5 lc rgb "#3465a4"
set style line 2 lt 1 lw 1.5 lc rgb "#f57900"
set xlabel "glutamate (mM)" offset 0,0.7;
fac=1.000;
set xrange [0.:0.06]
set yrange [-0.1:*]
set xtics nomirror 0.02;
set ytics nomirror;
set border 3;
eflow=3; sqo=eflow+1; fsc2=sqo+1; qsc2=fsc2+1; fsq1=qsc2+1; fum=fsq1+1; bypas=fum+1; qh=bypas+1; psi=qh+1; fc1=psi+1; rosc2=fc1+1; rosc3=rosc2+1; fc2=rosc3+1;
ca=fc2+1;
Glu_m=ca+1; OAAm=Glu_m+1; psio=OAAm+1; Na_i=psio+1; OAAc=Na_i+1; Glu_o=OAAc+1; Glu_i=Glu_o+1; ATP=Glu_i+1;
set key vert at graph 0.7,0.7 left Right samplen 1;
#set ylabel "ROS in C3, nM" offset 2;
#plot '00000' using ($1/fac):rosc3 w l ls 2 not;
#set ytics 0.2;
#set label 1 "A" at graph 0.05,1.0 font "arial,12"
#set ylabel "SQ at Qo in CIII (relative)" offset 2;
#plot '00000' using ($1/fac):sqo w l ls 2 not, 'a00000' using ($1/fac):sqo w l ls 1 not
set yrange [0:*]
#set ytics 20;
#set label 1 "B"
set label 1 "A" at graph 0.05,1.0 font "arial,12"
set label 2 "functional states" at graph 0.25,0.5 rotate by 25 
set label 3 "signaling states" at graph 0.25,0.95 
set label 4 "unstable states" at graph 0.25,0.87  rotate by -7
set ylabel "QH_2/(pool size)" offset 2;
plot '00000' using ($1/fac):qh w p ls 2 t"IVP", "mex/curveT.txt" u 15:($1/2.4) w l ls 1 t"MatContL", "../mark/noPTP/LP" u 7:($1/2.4) w p ps 1 pt 7 t"LP", "../mark/noPTP/hopf" u 7:($1/2.4) w p ps 1 lt 7 t"Hopf"
#set ylabel "ROS in C2, nM" offset 2;
#plot '00000' using ($1/fac):rosc2 w l ls 2 not;
#set yrange [-0.1:*]
#set ytics 0.2;
#set label 1 "C"
#set ylabel "ATP, mM" offset 2;
#plot '00000' using ($1/fac):ATP w l ls 2 not, 'a00000' using ($1/fac):ATP w l ls 1 not
set yrange [-5:*]
#set ytics 100;
set label 1 "B"
set label 2 at graph 0.3,1.0 rotate by -20
set label 3 "signaling" at graph 0.8,0.07 
set label 4 "unstable" at graph 0.57,0.07 rotate by 0
#set label 1 "D"
set ylabel "Ψ(var2, mV)" offset 2;
plot '00000' using ($1/fac):psi  w p ls 2 t"IVP", "mex/curveT.txt" u 15:($2*118.) w l ls 1 t"MatContL", "../mark/noPTP/LP" u 7:($2*118.) w p ps 1 pt 7 t"LP", "../mark/noPTP/hopf" u 7:($2*118.) w p ps 1 lt 7 t"Hopf"#, '00000' using ($1/fac):psi w l ls 1 not;
#, '00000' using ($1/fac):eflow w l ls 1 not;

