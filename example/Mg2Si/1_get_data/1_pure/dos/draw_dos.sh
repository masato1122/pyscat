FILE1="dos.txt"
FILE2="../../2_impurity/dos/dos.txt"

lw=3.0; ps=0.5;
gnuplot<<EOF
set size ratio 0.7 0.8

set xlabel "Frequency (THz)"
set xrange [*:*]
set xtics 2; set mxtics 2

set ylabel "DOS (a.u.)"
set yrange [*:*]
set ytics 1; set mytics 2

set out "out.eps"
plot \
    "$FILE1" u 1:((\$2)/2)  w lp ls 2 lw $lw ps $ps dt (0,0) title "pure",\
    "$FILE2" u 1:((\$2)/96) w lp ls 3 lw $lw ps $ps dt (0,0) title "impurity"
EOF

#EFIG=fig.eps
#ps2eps -l -f -B out.eps; mv out.eps.eps $EFIG; epstopdf $EFIG
#echo $EFIG; rm out.eps $EFIG

PFIG=fig_dos.png
convert -trim -transparent white -density 400 out.eps $PFIG
echo $PFIG
