export petsc_inputfile="/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/runs/petscrc.dat"
export mpiexec_petsc="/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/build/external/petsc/arch-linux2-c-opt/bin/mpiexec"
export petsc_solver="/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/build/ch2d_explicit"
#export petsc_solver="/home/adegennaro/Projects/AEOLUS/cahnhilliard_2d/cpp/build/ch2d_split"

echo "PETSc inputfile =" $petsc_inputfile
echo "PETSc solver =" $petsc_solver

# Run solver
$mpiexec_petsc -np 4 $petsc_solver $petsc_inputfile
#$mpiexec_petsc -np 4 valgrind --tool=memcheck -q --num-callers=20 --log-file=valgrind.log.%p $petsc_solver $petsc_inputfile
#valgrind --tool=memcheck -q --num-callers=20 --log-file=valgrind.log.%p $petsc_solver $petsc_inputfile
#$petsc_solver $petsc_inputfile
