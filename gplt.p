set xlabel "Time (s)" offset 0,0.7;
set xrange [0 to *]
set style line 1 lt 1 lw 1.5 lc rgb "#3465a4"
set style line 2 lt 1 lw 3 lc rgb "#f57900"
set style line 3 lt 1 lw 2 lc rgb "#4e9a06"
set style line 4 lt 1 lw 2 lc rgb "#ef2929"
eflow=2; sqo=eflow+1; fsc2=sqo+1; qsc2=fsc2+1; fsq1=qsc2+1; fum=fsq1+1; suc=fum+1; qh=suc+1; psi=qh+1; fc1=psi+1;
OAAm=fc1+1; naspm=OAAm+1; nglu=naspm+1; ATP=nglu+1; fc2=psi+1; ca=fc2+1; rosc2=fc1+1; rosc3=rosc2+1;

#//plot '011110rbm' u ($1):($2) every 3::3::293 ps 0.25 lc rgb "#3465a4" not;
set xtics 50 nomirror;
set ytics nomirror;
set border 3;
set key hor at graph 0.2,1.15 left Right samplen 1;

#set t svg enhanced dashed font "arial,10" size 700,2.5*68;
#set output fno;
#set multiplot layout 1,3;
#set label 1 "A" at graph 0.1,1.0 font "arial,12"
#set ytics 0.02
#set ylabel "FMN SQ in CI " offset 2;
#plot fn1 using ($1):fsq1 w l ls 1 not, fn2 using ($1):fsq1 w l ls 2 not
#set label 1 "B"
#set yrange [-0.045 to *]
#set ytics 0.1
#set ylabel "FAD SQ in CII " offset 2;
#plot fn1 using ($1):fsc2 w l ls 1 t"low", fn2 using ($1):fsc2 w l ls 2 t"high";
#set label 1 "C"
#set ytics 0.2
#set ylabel "SQ at Qo in CIII " offset 2;
#plot fn1 using ($1):sqo w l ls 1 not, fn2 using ($1):sqo w l ls 2 not;
#--------------------------------------------------------------

set t pngcairo enhanced dashed font "arial,10" size 700,5*68;
set output fno;
set multiplot layout 2,3;
set yrange [-0.0 to *]
#set ytics auto

set label 1 "A" at graph 0.07,1.05 font "arial,12"
set ytics 0.2
set ylabel "QH_2 " offset 2;
plot fn1 using ($1):qh w l ls 1 not, fn2 using ($1):qh w l ls 2 not;

set label 1 "B" 
set ytics 5
set ylabel "succinate (nmol/mg)" offset 1.5;
plot fn1 using ($1):suc w l ls 1 not, fn2 using ($1):suc w l ls 2 not;

set label 1 "C" 
set ytics 50
set ylabel "ΔΨ, vO_2" offset 2;
plot fn1 using ($1):eflow w l ls 1 t"vO_2", fn2 using ($1):psi w l ls 4 t"ΔΨ"#, fn2 using ($1):eflow w l ls 2 not, fn1 using ($1):psi w l ls 3 t"ΔΨ";

set label 1 "D" 
set ytics 0.2
set ylabel "ATP" offset 2;
plot fn1 using ($1):ATP w l ls 1 not, fn2 using ($1):ATP w l ls 2 not;

set label 1 "E"  
set ytics 0.002
set ylabel "OAA mito (mM)" offset 2;
plot fn1 using ($1):OAAm w l ls 1 not, fn2 using ($1):OAAm w l ls 2 not;

set label 1 "F" 
set ytics 1
set yrange [-0.0 to 3.3]
set ylabel "Glu mito (mM)" offset 2;
plot fn1 using ($1):nglu w l ls 1 not, fn2 using ($1):nglu w l ls 2 not;


#set label 1 "G"# at graph 0.05,1.0 font "arial,12"
#set ylabel "Asp mito (mM)" offset 1;
#plot fn1 using ($1):naspm w l ls 2 not#"work"# fn2 using ($1):rosc2 w l ls 3 t"glu";

#set label 1 "H" 
#set ylabel "Fum mito (mM)" offset 2;
#plot fn1 using ($1):fum w l ls 2 not#, fn2 using ($1):Glu_i w l ls 3 not, "rot1" using ($1):Glu_i w l ls 1 not;
#set ylabel "qsc2" offset 2;
#plot fn1 using ($1):qsc2 w l ls 1 not, fn2 using ($1):qsc2 w l ls 2 not;
unset multiplot


