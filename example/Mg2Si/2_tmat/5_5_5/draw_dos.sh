FILE="dos.txt"
FILE11="../../1_get_data/1_pure/dos/dos.txt"
FILE12="../../1_get_data/2_impurity/dos/dos.txt"

lw=1.0; ps=0.5;
gnuplot<<EOF
set size ratio 0.7 0.8
set key at graph 0.05,0.98 top left
set xlabel ""
set xrange [*:*]
#set xtics 1; set mxtics 2

set ylabel ""
set yrange [*:*]
#set ytics 1; set mytics 2

set out "out.eps"
plot \
    "$FILE11" u 1:((\$2)/96)        w p ls 2 lw $lw ps $ps pt 5 dt (0,0) title "pure (normal)",\
    "$FILE12" u 1:(3*(\$2)/125/125) w p ls 4 lw $lw ps $ps pt 9 dt (0,0) title "imp (normal)",\
    "$FILE"   u 1:2                 w p ls 1 lw $lw ps $ps pt 4 dt (0,0) title "pure",\
    "$FILE"   u 1:3                 w p ls 3 lw $lw ps $ps pt 8 dt (0,0) title "imp"
EOF

#EFIG=fig.eps
#ps2eps -l -f -B out.eps; mv out.eps.eps $EFIG; epstopdf $EFIG
#echo $EFIG; rm out.eps $EFIG

PFIG=fig_dos.png
convert -trim -transparent white -density 400 out.eps $PFIG
echo $PFIG
