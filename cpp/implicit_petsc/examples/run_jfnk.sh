export petsc_inputfile="petscrc.dat"

echo "PETSc inputfile =" $petsc_inputfile

# Run solver
./jfnk_2d $petsc_inputfile
#mpiexec -n 1 ./jfnk_2d $petsc_inputfile
