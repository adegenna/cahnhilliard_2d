export petsc_inputfile="petscrc.dat"

echo "PETSc inputfile =" $petsc_inputfile

# Run solver
./jfnk_2d $petsc_inputfile
