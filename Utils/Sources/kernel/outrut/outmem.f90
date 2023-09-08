!-----------------------------------------------------------------------
!> @addtogroup Memory
!> @{
!> @file    outmem.f90
!> @author  houzeaux
!> @date    2018-12-30
!> @brief   Output memory
!> @details Output information on memory required
!> @} 
!-----------------------------------------------------------------------

subroutine outmem()

  use def_master
  use def_domain
  use def_solver
  use def_inpout
  use def_kermod
  use mod_parall
  use mod_communications
  use mod_outfor, only : outfor
  use mod_memory, only : kfl_varcount
  use mod_memory, only : lun_varcount
  use mod_memory, only : memory_output_variable_counter
  implicit none
  real(rp)     :: rgiga,rmega,rkilo,rbyte
  integer(8)   :: imodu,iserv
  character(6) :: lbyte
  integer(ip)  :: ipass,number_passes

  real(rp)     :: r_tomax,r_memor_dom,r_memor_sol,r_memor_els
  real(rp)     :: r_mem_modul(mmodu)
  real(rp)     :: r_mem_servi(mserv)

  if( npart > 1 ) then
     number_passes = 2
  else
     number_passes = 1
  end if
  !
  ! First pass is for Master's max memory
  ! Second pass computes the max memory over the slaves
  !
  do ipass = 1,number_passes

     ioutp(50) = ipass     
     !
     ! Main memory: domain+master+solver
     !
     r_memor_dom = real(memor_dom(2),rp) 
     r_memor_sol = real(memma(2) + memdi(2) + memit(2),rp)     
     r_memor_els = relse(3)
     r_tomax     = r_memor_dom + r_memor_sol + r_memor_els
     do imodu = 1,mmodu
        if( kfl_modul(imodu) /= 0 ) then
           r_mem_modul(imodu) = real(mem_modul(2,imodu),rp)
           r_tomax = r_tomax + r_mem_modul(imodu)
        end if
     end do
     do iserv = 1,mserv
        if( kfl_servi(iserv) /= 0 ) then
           if( iserv == ID_PARALL ) then
             r_mem_servi(iserv) = real(mem_servi(2,iserv),rp) + real(par_memor(2),rp)
           else
              r_mem_servi(iserv) = real(mem_servi(2,iserv),rp)
           end if
           r_tomax = r_tomax + r_mem_servi(iserv)
        end if
     end do
     !
     ! Max values
     !
     if( ipass == 2 ) then 
        call PAR_MAX(r_tomax     , 'IN MY CODE')
        call PAR_MAX(r_memor_dom , 'IN MY CODE')
        call PAR_MAX(r_memor_sol , 'IN MY CODE')
        call PAR_MAX(r_memor_els , 'IN MY CODE')
        do imodu = 1,mmodu
           if( kfl_modul(imodu) /= 0 ) then
              call PAR_MAX(r_mem_modul(imodu),'IN MY CODE')
           end if
        end do
        do iserv = 1,mserv
           if( kfl_servi(iserv) /= 0 ) then
              call PAR_MAX(r_mem_servi(iserv),'IN MY CODE')
           end if
        end do
     end if
     !
     ! Gbutes, Mbytes or bytes?
     !
     rgiga = 1024_rp*1024_rp*1024_rp
     rmega = 1024_rp*1024_rp
     rkilo = 1024_rp     
     if( r_tomax >= rgiga ) then
        rbyte = rgiga
        lbyte = 'Gbytes'
     else if( r_tomax >= rmega ) then 
        rbyte = rmega
        lbyte = 'Mbytes' 
     else if( r_tomax >= rkilo ) then 
        rbyte = rkilo
        lbyte = 'kbytes'          
     else  
        rbyte = 1.0_rp
        lbyte = ' bytes'     
     end if

     routp(1) = r_tomax     / rbyte
     coutp(1) = lbyte
     routp(2) = r_memor_dom / rbyte
     coutp(2) = lbyte
     routp(3) = r_memor_els / rbyte
     coutp(3) = lbyte

     if( INOTSLAVE ) call outfor(21_ip,lun_outpu,' ')
     !
     ! Memory depending on the module
     !
     do imodu = 1,mmodu
        if( kfl_modul(imodu) /= 0 ) then
           coutp(1) = trim(namod(imodu))
           routp(1) = r_mem_modul(imodu) / rbyte
           coutp(2) = lbyte
           if( INOTSLAVE ) call outfor(22_ip,lun_outpu,' ')
        end if
     end do
     !
     ! Memory depending on the service
     !
     do iserv = 1,mserv
        if( kfl_servi(iserv) /= 0 ) then
           coutp(1) = trim(naser(iserv))
           routp(1) = r_mem_servi(iserv) / rbyte
           coutp(2) = lbyte
           if( INOTSLAVE ) call outfor(23_ip,lun_outpu,' ')
        end if
     end do
     !
     ! Solver and maximum memory
     !      
     routp(1) = r_memor_sol / rbyte
     coutp(1) = lbyte
     if( INOTSLAVE ) call outfor(24_ip,lun_outpu,' ')

  end do

  !----------------------------------------------------------------------
  !
  ! Variable memory counter
  !
  !----------------------------------------------------------------------

  if( kfl_varcount == 1 ) then
     call memory_output_variable_counter(lun_varcount,OUTPUT_FORMAT="('VARIABLE= ',a,' MEMORY= ',e13.6)")
  end if

end subroutine outmem
