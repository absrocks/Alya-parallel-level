subroutine ale_smodef()
  !-----------------------------------------------------------------------
  !****f* domain/ale_smodef
  ! NAME
  !    domain
  ! DESCRIPTION
  !    This routines 
  !    KFL_FIXNO_ALE(:,:) = -1 ... FMALE free
  !                       =  0 ... ALE free
  !                       =  1 ... ALE fixed
  !                       =  3 ... FMALE fixed
  ! USED BY
  !    Turnon 
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_elmtyp
  use def_domain
  use def_alefor
  use mod_memchk
  use def_kermod,         only : kfl_adj_prob
  use mod_messages, only : livinf

  implicit none
  integer(ip)    :: idime,ipoin
  character(150) :: messa

  if( kfl_defor_ale /= 0 .or. kfl_smoot_ale /= 0 .or. kfl_smobo_ale /= 0 ) then

     !-------------------------------------------------------------------
     !
     ! Copy new coordinates
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then
        do ipoin = 1,npoin
           do idime = 1,ndime
              coord_ale(idime,ipoin,1) = coord(idime,ipoin) 
           end do
        end do
     end if

     !-------------------------------------------------------------------
     !
     ! Boundary smoothing: Sharp edges have KFL_FIXNO_ALE(1,IPOIN) = 2
     !
     !-------------------------------------------------------------------

     if( kfl_smobo_ale /= 0 ) then
        call runend('CHECK ASSEMBLY, NO LONGER SUPPORTED')
        messa = intost(nsmob_ale) 
        messa = 'BOUNDARY SMOOTHING USING '//trim(messa)//' STEPS'
        call livinf(79_ip,trim(messa),modul)
        call ale_smobou()
     end if

     !-------------------------------------------------------------------
     !                   
     ! Mesh deformation
     !
     ! Update displacement and mesh velocity
     ! FMALE: 1. does not update coordinate
     !        2. Invert mesh velocity
     !  
     !-------------------------------------------------------------------
 
     if( kfl_defor_ale /= 0 ) then
        messa = intost(ndefo_ale) 
        messa = 'MESH DEFORMATION USING '//trim(messa)//' LOADING STEPS'
        call livinf(79_ip,trim(messa),modul)
        call ale_deform()
        call ale_updunk(3_ip)
     end if

!!$  if( INOTMASTER ) then
!!$     do iimbo = 1,nrbod
!!$        do kpoin = 1,rbbou(iimbo) % npoib
!!$           ipoin = rbbou(iimbo) % lninv(kpoin)
!!$           do idime = 1,ndime    
!!$              bvess_ale(idime,ipoin)     = rbbou(iimbo) % cooib(idime,kpoin) - coord(idime,ipoin)
!!$              kfl_fixno_ale(idime,ipoin) = 1
!!$           end do
!!$        end do
!!$     end do
!!$  end if 

      !if( INOTMASTER ) then
      !   do iperi = 1,nperi
      !       ipoin = lperi(1,iperi)
      !       jpoin = lperi(2,iperi)
      !       if ( ipoin > 0 .and. jpoin > 0 ) then
      !          do idime = 1,ndime
      !             if( abs(dispm(idime,ipoin,1)-dispm(idime,jpoin,1)) > 1.0e-12_rp ) then
      !                write(*,*) 'MERDE 777=',abs(dispm(idime,ipoin,1)-dispm(idime,jpoin,1))
      !             end if
      !          end do
      !       end if
      !    end do
      ! end if
      ! call par_check_continuity(ndime,dispm)

     !-------------------------------------------------------------------
     !                                                    
     ! Mesh smoothing
     !
     ! Update displacement and mesh velocity
     ! FMALE: 1. does not update coordinate
     !        2. Invert mesh velocity
     !                                                                   
     !-------------------------------------------------------------------

     if( kfl_smoot_ale /= 0 ) then
        messa = intost(nsmoo_ale) 
        messa = 'MESH SMOOTHING USING '//trim(messa)//' STEP'
        call livinf(79_ip,trim(messa),modul)
        call ale_smooth()
        call ale_updunk(3_ip)
     end if

     !-------------------------------------------------------------------
     ! 
     ! Update displacement and mesh velocity
     ! FMALE: 1. does not update coordinate
     !        2. Invert mesh velocity
     !
     !-------------------------------------------------------------------
     !call ale_updunk(3_ip)

  end if
  
  if (kfl_adj_prob == 1) then
    
     !-------------------------------------------------------------------
     ! 
     ! Update mesh sensitivities 
     ! 
     !-------------------------------------------------------------------

     call ale_mshsen
     return
  endif
  

end subroutine ale_smodef
