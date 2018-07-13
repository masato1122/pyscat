FILE1="dos.txt"
FILE2="../../1_get_data/1_pure/dos/dos_555.txt"
FILE3="../../1_get_data/1_pure/dos_green/dos.txt"

lw=3.0; ps=0.5;
gnuplot<<EOF
set size ratio 0.7 0.8
set key at graph 0.05,0.95 left top

set xlabel ""
set xrange [*:*]
#set xtics 1; set mxtics 2

set ylabel ""
set yrange [*:*]
#set ytics 1; set mytics 2

set out "out.eps"
plot \
    "$FILE1" u 1:2 w p ls 2 lw $lw ps $ps title "DOS(green1)",\
    "$FILE3" u 1:((\$2)/6) w p ls 3 lw $lw ps $ps title "DOS(green2)"
EOF

#"$FILE2" u 1:((\$2)/96) w p ls 4 lw $lw ps $ps title "DOS(normal)"

#EFIG=fig.eps
#ps2eps -l -f -B out.eps; mv out.eps.eps $EFIG; epstopdf $EFIG
#echo $EFIG; rm out.eps $EFIG

PFIG=fig_dos.png
convert -trim -transparent white -density 400 out.eps $PFIG
echo $PFIG
