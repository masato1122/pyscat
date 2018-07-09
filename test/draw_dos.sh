#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)

FILE="nacl_dos.txt"

lw=3.0; ps=0.5;
gnuplot<<EOF
set xlabel "Frequency (THz)"
set xrange [*:*]
#set xtics 1; set mxtics 2

set ylabel "DOS (au)"
set yrange [*:*]
#set ytics 1; set mytics 2

set out "out.eps"
plot \
    "$FILE" u 1:2 w p ls 2 lw $lw ps $ps title ""
EOF

#EFIG=fig.eps
#ps2eps -l -f -B out.eps; mv out.eps.eps $EFIG; epstopdf $EFIG
#echo $EFIG; rm out.eps $EFIG

PFIG=NaCl_dos.png
convert -trim -transparent white -density 400 out.eps $PFIG
echo $PFIG
