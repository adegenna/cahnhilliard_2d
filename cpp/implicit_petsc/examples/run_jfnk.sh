export petsc_inputfile="petscrc.dat"

echo "PETSc inputfile =" $petsc_inputfile

# Run solver
mpiexec -np 8 ./jfnk_2d $petsc_inputfile
