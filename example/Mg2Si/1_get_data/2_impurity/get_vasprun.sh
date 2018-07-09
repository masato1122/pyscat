
DIR0=/work/k0247/k024705/Impurity/example/Mg2Si/2_impurity
for i in `seq 1 20`; do
    scpsekirei_get $DIR0/disp-${i}/vasprun.xml
    mv ./vasprun.xml ./disp-${i}/vasprun.xml
done

