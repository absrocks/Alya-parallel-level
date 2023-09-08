subroutine nsa_outset()
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_outset
  ! NAME 
  !    nsa_outset
  ! DESCRIPTION
  !    Compute and write results on sets:
  !    - Element sets:
  !      1.    VELOC: mean velocity =  sqrt{ (int_W u^2 dw)/(2*meas(W)) }
  !      2.    VORTI: mean vorticity = sqrt{ (int_W w^2 dw)/(2*meas(W)) }
  !                                    with w=dv/dx-du/dy
  !      3.    VORTI: mean vorticity = int_W rho*u^2/2 dw
  !
  !    - Boundary sets:
  !      1.    MEANP: mean pressure  =  int_W p/meas(W)
  !      2.    MASS:  mass           =  int_W rho*u.n 
  !      3-8.  FORCE: force          =  int_W [-pI+2 mu E(u)].n 
  !      9-14. TORQU: torque         =  int_W r x ( [-pI+2 mu E(u)].n )
  !      15.   MEANY: mean y+        =  int_W y+/meas(W)
  !      Force and torque are exerted by the solid on the fluid.
  !
  !   - Node sets:
  !      1-3.  VELOC: velocity components
  !      4.    PRESS: pressure
  !      5.    YPLUS: y+
  ! USES
  ! USED BY
  !    nsa_output
  !***
  !----------------------------------------------------------------------- 
  use def_parame
  use def_master
  use def_domain
  use mod_iofile

  use def_nastal

  implicit none
  integer(ip) :: ieset,ibset,inset,idime,dummi
  real(rp)    :: resub,ubulk,ubdif

  !----------------------------------------------------------------------
  !
  ! Element sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setse) > 0 ) then

     if( INOTMASTER ) then
        do ieset = 1,neset
           call nsa_elmset(lesec(ieset),ieset)
        end do
     end if
     call posdef(21_ip,dummi)
     !
     ! Normalize results of necessary
     !
     if( INOTSLAVE ) then
        do ieset = 1,neset
           if(postp(1) % npp_setse(1)/=0.and.veset(postp(1) % nvaes+1,ieset)>0.0_rp) then
              veset(1,ieset)=sqrt(veset(1,ieset)/(2.0_rp*veset(postp(1) % nvaes+1,ieset)))
           end if
           if(postp(1) % npp_setse(2)/=0.and.veset(postp(1) % nvaes+1,ieset)>0.0_rp) then
              veset(2,ieset)=sqrt(veset(2,ieset)/(2.0_rp*veset(postp(1) % nvaes+1,ieset)))
           end if
        end do
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Boundary sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsb) > 0 ) then

     if( INOTMASTER ) then
        do ibset = 1,nbset 
           call nsa_bouset(lbsec(ibset),ibset)
        end do
     end if

     call posdef(22_ip,dummi) !! It inlcudes an MPI_allreduce

     !
     ! Normalize results of necessary
     !
     do ibset = 1,nbset
        if(postp(1) % npp_setsb(1)/=0.and.vbset(postp(1) % nvabs+1,ibset)>0.0_rp) then
           vbset( 3,ibset) = - vbset( 1,ibset) / vbset( 2,ibset)
        end if
     end do

     !
     ! Mass flow rate control
     !
     if( kfl_mfrco_nsa > 0 ) then

       if (mod(ittim,10_ip) == 0) then

          if ( INOTMASTER ) then
    
             if (ittim == 1) ubpre_nsa = vbset( 3,mfrse_nsa)   !! DMM This condition must be changed for restart

             resub = mfrub_nsa - vbset( 3, mfrse_nsa)
             ubdif = vbset( 3, mfrse_nsa) - ubpre_nsa
                  
             if ( (resub / mfrub_nsa) > 0.03_rp .and. ubdif < 0) then
                do idime=1,ndime
                   gravm_nsa(idime,ndime+1) = mfrgr_nsa*gravm_nsa(idime,ndime+1)
                   gravm_nsa(ndime+2,idime) = mfrgr_nsa*gravm_nsa(ndime+2,idime)
                end do
             elseif ( (resub / mfrub_nsa) < -0.03_rp .and. ubdif > 0) then
                do idime=1,ndime
                   gravm_nsa(idime,ndime+1) = (2.0_rp - mfrgr_nsa) * gravm_nsa(idime,ndime+1)
                   gravm_nsa(ndime+2,idime) = (2.0_rp - mfrgr_nsa) * gravm_nsa(ndime+2,idime)
                end do
             endif
         
             write(11111,*) ittim,vbset( 3, mfrse_nsa),gravm_nsa(1,ndime+1)

             ubpre_nsa = vbset( 3, mfrse_nsa)

          endif
       endif

     endif

  end if

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsn) > 0 ) then

     if( INOTMASTER ) then
        do inset = 1,nnset
           if(lnsec(inset)/=0) then
              do idime=1,ndime
                 if(postp(1) % npp_setsn(idime)/=0)&
                      vnset(idime,inset)=veloc(idime,lnsec(inset),1)
              end do
              if(ndime==2) vnset(3,inset)=0.0_rp
              if(postp(1) % npp_setsn(4)/=0) vnset(4,inset)=press(lnsec(inset),1)
           end if
        end do
     end if
     call posdef(23_ip,dummi)

  end if

end subroutine nsa_outset
