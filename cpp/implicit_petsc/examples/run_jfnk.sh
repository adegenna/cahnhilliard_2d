export petsc_inputfile="petscrc.dat"
export mpiexec_petsc="/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/build/external/petsc/arch-linux2-c-opt/bin/mpiexec"

echo "PETSc inputfile =" $petsc_inputfile

# Run solver
$mpiexec_petsc -np 1 ./jfnk_2d $petsc_inputfile
