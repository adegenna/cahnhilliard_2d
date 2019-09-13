export petsc_inputfile="/home/adegennaro/cahnhilliard_2d/cpp/implicit_petsc/examples/petscrc.dat"
export mpiexec_petsc="/home/adegennaro/cahnhilliard_2d/cpp/build/external/petsc/arch-linux2-c-opt/bin/mpiexec"

echo "PETSc inputfile =" $petsc_inputfile

# Run solver
$mpiexec_petsc ./jfnk_2d $petsc_inputfile
