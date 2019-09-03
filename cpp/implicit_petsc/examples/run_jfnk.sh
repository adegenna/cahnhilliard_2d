export petsc_inputfile="petscrc.dat"

echo "PETSc inputfile =" $petsc_inputfile

# Run solver
#./jfnk_2d $petsc_inputfile
mpiexec -np 4 ./jfnk_2d $petsc_inputfile
