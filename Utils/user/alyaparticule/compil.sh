#
#  compil all with mpif90
#
mpif90 alya-deposition.f90 -o alya-deposition.x
mpif90 alya-particule.f90 -o alya-particule.x
mpif90 denis_depo.f90 -o denis_depo.x
mpif90 depo-initial.f90 -o depo-initial.x
mpif90 slice-particule.f90 -o slice-particule.x
mpif90 valid_nasal_deposition.f90  -o valid_nasal_deposition.x
mpif90 valid_pipe_deposition.f90  -o valid_pipe_deposition.x
mpif90 slice_nasal1_deposition.f90  -o slice_nasal1_deposition.x
cp alya-particule.x ../.
cp alya-deposition.x ../.
cp valid_nasal_deposition.x ../.
cp valid_pipe_deposition.x ../.
cp slice_nasal1_deposition.x ../.
cp denis_depo.x ../.

