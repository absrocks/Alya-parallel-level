subroutine outhel
  !------------------------------------------------------------------------
  !****f* master/outhel
  ! NAME
  !    outhel
  ! DESCRIPTION
  !    This routine output a help message
  ! OUTPUT
  !   
  ! USED BY
  !    openfi
  !***
  !------------------------------------------------------------------------

  use      def_master
  implicit none

  write(6,1)
  stop

1 format(&
       & '   ',/,&
       & '   ALYA USAGE:',/,&
       & '   ',/,&
       & '   Alya.x [-h -c] name ',/,&
       & '   ',/,&
       & '   -h --help  ......... Displays this help',/,&
       & '   -c --check ......... Check data file',/,&
       & '   -f --file name ..... Problem name',/,&
       & '   ',/,&
       & '   ',/,&
       & '   Runs Alya for problem name. ',/,&
       & '   ',/,&
       & '   The following I/O files are located/created in current directory (mod is any activated module extension)',/,&
       & '   * means optional:',/,&
       & '   ',/,&
       & '   (I)    name.dat:                     run data',/,&
       & '   (I)    name.ker.dat:                 kernel data',/,&
       & '   (I)    name.dom.dat:                 mesh data',/,&
       & '   (I*)   name.cou.dat:                 coupling data',/,&
       & '   (I)    name.mod.dat:                 module data',/,&
       & '   ',/,&
       & '   (O)    name.log:                     run log',/,&      
       & '   (O)    name.ker.log:                 kernel log',/,&      
       & '   (O*)   name.mem:                     memory',/,&
       & '   (O*)   name.liv:                     live info',/,&
       & '   (O)    name-partition.par.post.msh   partition mesh in GiD format',/,&
       & '   (O)    name-partition.par.post.res   partition results in GiD format',/,&
       & '   (O)    name-VAR.post.alyabin:        postprocess file of variable VAR',/,&
       & '   (O)    name-VAR.mod.sol              solver information for variable VAR',/,&
       & '   (O)    name-VAR.mod.cso              solver convergence for variable VAR',/,&
       & '   (O)    name.mod.cvg:                 module convergence',/,&      
       & '   (O)    name.mod.rst:                 module restart',/,&  
       & '   ')

end subroutine outhel
