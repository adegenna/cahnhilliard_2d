### MAKE INITIAL DATA FILES USING PYTHON SCRIPT
python make_input_data_files_test1.py


### RUN THE ACTUAL CH SOLVER
export build_dir="../../build/"

export mpiexec_petsc=${build_dir}"external/petsc/arch-linux2-c-opt/bin/mpiexec"
export petsc_inputfile="./petscrc.dat"
#export petsc_solver=${build_dir}"ch2d_split"
export petsc_solver=${build_dir}"ch2d_explicit"

echo "PETSc inputfile =" $petsc_inputfile
echo "PETSc solver =" $petsc_solver

# Run solver
$mpiexec_petsc -np 4 $petsc_solver $petsc_inputfile
#$mpiexec_petsc -np 4 valgrind --tool=memcheck -q --num-callers=20 --log-file=valgrind.log.%p $petsc_solver $petsc_inputfile
#valgrind --tool=memcheck -q --num-callers=20 --log-file=valgrind.log.%p $petsc_solver $petsc_inputfile
#$petsc_solver $petsc_inputfile
