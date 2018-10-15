from ase.io import (read, write)

#type = "pure"
type = "imp"
pfile = "SPOSCAR_" + type
pxyz = type + ".xyz"

atoms = read(pfile, index=0, format="vasp")
write(pxyz, images=atoms, vec_cell=True)
write('POSCAR1', images=atoms, format='vasp',
        direct=True, vasp5=True, sort=True)


