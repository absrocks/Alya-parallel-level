subroutine qua_outinf()
!-----------------------------------------------------------------------
  !****f* Quanty/qua_outinf
  ! NAME 
  !    qua_outinf
  ! DESCRIPTION
  !    This routine writes on the quanty files
  ! USES
  !    qua_inisol
  !    outfor
  !    outsol
  ! USED BY
  !    qua_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_quanty
  implicit none
  
  character(60) :: equat
  integer(ip)   :: ierhs,imate
  character(2)  :: wcpcv

  if(kfl_paral<=0) then
     !
     ! Write information in Result file
     !
     if(kfl_rstar/=2) then
        equat=''
     end if

  end if
  !
  ! Formats
  !
110 format(/,&
       & 10x,a,//,&
       & 10x,'Physical properties: ',i2,' material(s)',/,&
       & 10x,'------------------------------------------------------------ ')
111 format(&
       & 10x,a,10(e12.6,1x))
112 format(&
       & 10x,a,e12.6,a)

end subroutine qua_outinf

