subroutine outres()
  !------------------------------------------------------------------------
  !****f* output/outres
  ! NAME 
  !    outres
  ! DESCRIPTION
  !    This routine composes mesh and results file names
  ! OUTPUT
  ! USES
  ! USED BY
  !    openfi
  !***
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_postpr
  implicit none

  if(kfl_outfo==1) then
     !
     ! GiD
     !
     fil_outpu_dom = trim(fil_postp)//'.post.msh'
     fil_postp     = trim(fil_postp)//'.post.res'

  else if(kfl_outfo==2) then
     !
     ! Femview
     !
     fil_outpu_dom = trim(fil_postp)//'.fem'
     fil_postp     = trim(fil_postp)//'.fem'

  else if(kfl_outfo>9.and.kfl_outfo<21) then
     !
     ! Ensight/Visit
     !
     fil_outpu_dom = trim(fil_postp)//'.ensi.geo'
     fil_pos00     = trim(fil_postp)//'.ensi.case'
     fil_postp     = trim(fil_postp)//'.ensi'

  else if(kfl_outfo==30) then
     !
     ! VU
     !
     fil_outpu_dom = trim(fil_postp)//'.msh.vu'
     fil_postp     = trim(fil_postp)//'.res'

  end if

end subroutine outres
