cp ./coord/POSCAR-Mg2Si-Sn ./POSCAR

FDIR=./FILES

TYPE=imp

# -- prepare directories for force calculations
phonopy -d --dim="1 1 1"
for i in `seq 1 50`; do
    NUM=`printf %03d ${i}`
    PFILE=POSCAR-${NUM}
    if [ ! -e $PFILE ]; then
        continue
    fi
    DIR=./disp-$i
    mkdir $DIR
    mv $PFILE $DIR/POSCAR
    for f in POTCAR KPOINTS INCAR; do
        cp $FDIR/$f $DIR
    done

    #---- for me
    DIR1=../UTILS
    for ss in dalton sekirei; do
        cp $DIR1/${ss}.sh ./$DIR
        sed -i -e "s/test/${TYPE}${NUM}/g" ./$DIR/${ss}.sh
    done
done
