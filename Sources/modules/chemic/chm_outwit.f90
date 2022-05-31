subroutine chm_outwit()
  !------------------------------------------------------------------------
  !****f* chemic/chm_outwit
  ! NAME 
  !    chm_outwit
  ! DESCRIPTION
  !    Output chemic species on witness points
  ! USES
  ! USED BY
  !    chm_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_chemic 
  use mod_memory,      only     : memory_alloca, memory_deallo
  use mod_interp_tab,  only     : fw_scale_cont_var 
  use def_chemic,      only     : table_fw

  implicit none
  integer(ip)       :: iwitn,ielem,inode,pnode,pelty,ipoin,iclas,ivawi,dummi
  
  real(rp) , pointer        :: control(:)   ! input of table lookup function 
  real(rp) , pointer        :: scale_control(:)
  real(rp) , pointer        :: lim_control(:,:)

  if( nwitn > 0 .and. maxval(postp(1) % npp_witne) > 0 ) then
     !
     ! Results on witness points
     !
     witne => postp(1) % witne
 
     if( INOTMASTER ) then

       !
       ! initialization for SCONC
       !
       if( postp(1) % npp_witne(8+8) == 1 ) then
           nullify(control)
           nullify(scale_control)
           nullify(lim_control)
           call memory_alloca(mem_modul(1:2,modul),'control',      'chm_outvar',control,      5_ip)
           call memory_alloca(mem_modul(1:2,modul),'scale_control','chm_outvar',scale_control,5_ip)
           call memory_alloca(mem_modul(1:2,modul),'lim_control',  'chm_outvar',lim_control,  5_ip, 2_ip)
       endif


       do iwitn = 1,nwitn
          ielem = lewit(iwitn)
          if( ielem > 0 ) then
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             do ivawi = 1,postp(1) % nvawi
                witne(ivawi,iwitn) = 0.0_rp
             end do 
      
             if( postp(1) % npp_witne(1) == 1 ) then
                do iclas = 1,min(8_ip,nclas_chm)
                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      witne(iclas,iwitn) = witne(iclas,iwitn) + shwit(inode,iwitn) * conce(ipoin,iclas,1) * postp(1) % witne_dt(1)
                   end do
                end do
             end if             

             if( postp(1) % npp_witne(8+1) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+1,iwitn) = witne(8+1,iwitn) + shwit(inode,iwitn) * xZr_chm(ipoin) * postp(1) % witne_dt(9)
                end do
             end if

             if( postp(1) % npp_witne(8+2) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+2,iwitn) = witne(8+2,iwitn) + shwit(inode,iwitn) * xZs_chm(ipoin) * postp(1) % witne_dt(10)
                end do
             end if

             if( postp(1) % npp_witne(8+3) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+3,iwitn) = witne(8+3,iwitn) + shwit(inode,iwitn) * xYr_chm(ipoin) * postp(1) % witne_dt(11) 
                end do
             end if

             if( postp(1) % npp_witne(8+4) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+4,iwitn) = witne(8+4,iwitn) + shwit(inode,iwitn) * xYs_chm(ipoin) * postp(1) % witne_dt(12)
                end do
             end if

             if( postp(1) % npp_witne(8+5) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+5,iwitn) = witne(8+5,iwitn) + shwit(inode,iwitn) * d32_chm(ipoin) * postp(1) % witne_dt(13) 
                end do
             end if

             if( postp(1) % npp_witne(8+6) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+6,iwitn) = witne(8+6,iwitn) + shwit(inode,iwitn) * Sigma_chm(ipoin)  * postp(1) % witne_dt(14)
                end do
            end if

            if( postp(1) % npp_witne(8+7) == 1 ) then

                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   witne(8+7,iwitn) = witne(8+7,iwitn) + shwit(inode,iwitn) * Sigm0_chm(ipoin) * postp(1) % witne_dt(15)
                end do
             end if
             
             
             if( postp(1) % npp_witne(8+8) == 1 ) then
                !
                ! Scaled control variables
                ! Assume CONCE is output as well, so the control variables are
                ! interpolated already except therm
                !
                control = 0.0_rp
                do iclas = 1, table_fw % main_table % ndim
                   select case (table_fw % main_table % coords(iclas) % name)
                   case ('CMEAN','C    ')
                       control(iclas) = witne(1,iwitn)
                   case ('CVAR ')
                       control(iclas) = witne(2,iwitn)
                   case ('CHIST')
                       control(iclas) = witne(8+1,iwitn) +  witne(8+2,iwitn)
                   case ('ZMEAN','Z    ')
                       control(iclas) = witne(3,iwitn)
                   case ('ZVAR ')
                       control(iclas) = witne(4,iwitn)
                   case ('IMEAN','I    ')
                       do inode = 1,pnode
                          ipoin = lnods(inode,ielem)
                          control(iclas) = control(iclas) + shwit(inode,iwitn) * therm(ipoin,1) * postp(1) % witne_dt(16)
                       enddo
                   end select
                enddo

                call fw_scale_cont_var( control, scale_control, lim_control, table_fw)

                do iclas = 1, table_fw % main_table % ndim
                    witne(8+7+iclas,iwitn) = scale_control(iclas) 
                end do
             end if             


          end if

       end do
       

       if( postp(1) % npp_witne(8+8) == 1 ) then
          call memory_deallo(mem_modul(1:2,modul),'control'      ,'chm_outvar',control)
          call memory_deallo(mem_modul(1:2,modul),'scale_control','chm_outvar',scale_control)
          call memory_deallo(mem_modul(1:2,modul),'lim_control'  ,'chm_outvar',lim_control)
       endif

     end if

  end if

end subroutine chm_outwit

