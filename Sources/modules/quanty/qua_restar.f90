subroutine qua_restar(itask)
  !------------------------------------------------------------------------
  !****f* Quanty/qua_restar
  ! NAME 
  !    qua_restar
  ! DESCRIPTION
  !    This routine reads the initial values from the restart file
  ! USES
  ! USED BY
  !    qua_turnon
  !***
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_quanty
  use      mod_iofile
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,icomp

  if(kfl_paral==-1) then

     select case (itask)

     case (1)
        if(kfl_rstar>=1) then
   !        call iofile(zero,lun_rstar_qua,fil_rstar_qua,'RESTART','old','unformatted') 
   !        icomp=min(3,ncomp_qua) 
   !        read(lun_rstar_qua) & 
   !             (phi(ipoin,icomp),ipoin=1,npoin)   
   !        call iofile(two,lun_rstar_qua,fil_rstar_qua,'RESTART') 
   !        call qua_updunk(six)
        end if
     case (2)
   !     if(kfl_preli==1.and.(&
   !          mod(ittim,nprit)==0.or.&
   !          cutim>=timef-zetem.or.&
   !          ittim>=mitim.or.&
   !          kfl_timei==0)) then
   !        call iofile(zero,lun_rstar_qua,fil_rstar_qua,'RESTART','unknown','unformatted')        
   !        write(lun_rstar_qua) &
   !             (tempe(ipoin,1),ipoin=1,npoin)
   !        call iofile(two,lun_rstar_qua,fil_rstar_qua,'RESTAR') 
   !     end if
     end select

  end if

end subroutine qua_restar
