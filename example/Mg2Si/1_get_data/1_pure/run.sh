#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)

for i in `seq 1 50`; do
    NUM=`printf %03d ${i}`
    DIR=$WORKDIR/disp-$i
    if [ ! -e $DIR ]; then
        continue
    fi
    cd $DIR
    qsub sekirei.sh
done
