!-----------------------------------------------------------------------
!> @addtogroup IO
!> @{
!> @file    mod_outvar.f90
!> @author  houzeaux
!> @date    2020-08-27
!> @brief   Output variables 
!> @details Output variables for time step jttim and time value dutim
!-----------------------------------------------------------------------

module mod_outvar

  use def_domain 
  use mod_memory
  use def_kintyp_basic,      only : ip,rp,lg
  use def_kintyp_mesh_basic, only : mesh_type_basic
  use def_parame,            only : zero,one
  use def_master,            only : ivapo,ITASK_ENDRUN,INOTMASTER,iasca,iavec,iar3p
  use def_master,            only : IMASTER,mitim,ittyp,intost,kfl_paral
  use def_master,            only : gescx,gesca,gevec,ger3p,gevex,postp
  use def_kermod,            only : witness_mesh
  use mod_postpr,            only : postpr
  use def_postpr,            only : kfl_ivari
  use mod_communications,    only : PAR_MAX
  use mod_arrays,            only : arrays
  use mod_interp_fe,         only : interp_fe
  use mod_interp_fe,         only : interp_fe_deallocate
  implicit none

  private

  public :: outvar
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-08-27
  !> @brief   Output
  !> @details Main output subroutine
  !> 
  !-----------------------------------------------------------------------


  subroutine outvar(jvari,jttim,dutim,wopos,MESH,MESH_ID)

    integer(ip),                      intent(in) :: jvari
    integer(ip),                      intent(in) :: jttim
    real(rp),                         intent(in) :: dutim
    character(5),                     intent(in) :: wopos(3)
    type(mesh_type_basic),  optional, intent(in) :: MESH
    integer(ip),            optional, intent(in) :: MESH_ID    
    integer(ip)                                  :: iosca,iovec,ior3p,itens,ipoin,ii
    integer(ip)                                  :: idime,ivari,kdime,dim_max,imesh
    character(5)                                 :: wauxi(3,6)
    character(2)                                 :: waux2(9)
    character(5)                                 :: waux3(3)
    character(20)                                :: wtype
    real(rp),               pointer              :: dusca(:)
    real(rp),               pointer              :: duvec(:,:)
    type(r3p),              pointer              :: dur3p(:)
    real(rp),               pointer              :: ggggg(:)
    real(rp),               pointer              :: vvvvv(:,:)
    real(rp),               pointer              :: gesca_loc(:)
    real(rp),               pointer              :: gevec_loc(:,:)
    logical(lg)                                  :: if_interpolated
    
    nullify(vvvvv)
    nullify(ggggg)
    nullify(dusca)
    nullify(duvec)
    nullify(dur3p)
    nullify(gesca_loc)
    nullify(gevec_loc)

    ivari = abs(jvari)
    ivapo = ivari
    if( mitim == 0 .and. ittyp == ITASK_ENDRUN ) goto 10

    iosca = 0
    iovec = 0
    ior3p = 0
    !
    ! A mesh is present
    !
    if( present(MESH_ID) ) then
       imesh = MESH_ID
    else if( present(MESH) ) then
       imesh = MESH % id
    else
       imesh = 0
    end if
    !
    ! Array should be interpolated 
    !
    if_interpolated = .false.
    if( imesh /= 0 ) then
       if( associated(witness_mesh(imesh) % inte % lelem) ) then
          if_interpolated = .true.
       end if
    end if

    if( wopos(2) == 'SCALA' .or. wopos(2) == 'SCALX' ) then
       !
       ! Postprocess scalar
       !
       iosca = 1

    else if ( wopos(2) == 'VECTO' .or. wopos(2) == 'VECTX' ) then
       !
       ! Postprocess vector
       !
       wauxi(2,1) = 'VECTO'
       wauxi(3,1) = wopos(3)    
       iovec = 1

    else if ( wopos(2) == 'R3P  ' .or. wopos(2) == 'R3PVE' ) then
       !
       ! Values at Gauss points
       !
       ior3p = 1

    else if ( wopos(2) == 'MULTI' ) then
       !
       ! Postprocess multidimensional array
       !
       iosca = 1

    else if ( wopos(2) == 'TENSO' ) then
       !
       ! Postprocess tensor: transform it into scalar
       !
       iovec = 1

    end if

    if ( wopos(2) == 'SCALA' ) then
       !
       ! SCALAR
       !
       if( kfl_ivari(2) > 0 ) then
          call posvox(wopos)
       end if
       if( kfl_ivari(1) > 0 ) then
          if( IMASTER ) then
             if( imesh == 0 ) then
                call postpr(dusca,wopos,jttim,dutim)
             else
                call postpr(dusca,wopos,jttim,dutim,MESH=witness_mesh(imesh) % mesh)
             end if
          else
             if( imesh == 0 ) then
                call postpr(gesca,wopos,jttim,dutim)
             else
                if( if_interpolated ) then
                   call interp_fe(gesca,witness_mesh(imesh) % inte,gesca_loc,witness_mesh(imesh) % mesh % npoin)
                   call postpr(gesca_loc,wopos,jttim,dutim,MESH=witness_mesh(imesh) % mesh)
                else
                   call postpr(gesca,wopos,jttim,dutim,MESH=witness_mesh(imesh) % mesh)
                end if
             end if
          end if
       end if

    else if ( wopos(2) == 'VECTO' ) then
       !
       ! VECTOR
       !
       if( kfl_ivari(2) > 0 ) then
          call posvox(wopos)
       end if
       if( kfl_ivari(1) > 0 ) then
          if( IMASTER ) then
             if( imesh == 0 ) then
                call postpr(duvec,wopos,jttim,dutim,ndime)
             else
                call postpr(duvec,wopos,jttim,dutim,ndime,MESH=witness_mesh(imesh) % mesh)              
             end if
          else
             if( imesh == 0 ) then
                call postpr(gevec,wopos,jttim,dutim,ndime)
             else
                if( if_interpolated ) then
                   call interp_fe(gevec,witness_mesh(imesh) % inte,gevec_loc,witness_mesh(imesh) % mesh % npoin)
                   call postpr(gevec_loc,wopos,jttim,dutim,MESH=witness_mesh(imesh) % mesh)
               else
                   call postpr(gevec,wopos,jttim,dutim,ndime,MESH=witness_mesh(imesh) % mesh)
                end if
             end if
          end if
       end if

    else if ( wopos(2) == 'R3P  ' .or. wopos(2) == 'R3PVE' ) then
       !
       ! R3P
       !
       if( kfl_ivari(1) > 0 ) then 
          if( imesh == 0 ) then
             if( IMASTER ) then
                call postpr(dur3p,wopos,jttim,dutim)
             else
                call postpr(ger3p,wopos,jttim,dutim)
             end if
          else
             call runend('OUTVAR: CANNOT POSTPROCESS R3P VARIABLES ON WTNESS MESHES')
          end if
       end if

    else if ( wopos(2) == 'MULTI' ) then
       !
       ! Multi-dimensional
       !
       if( kfl_ivari(1) > 0 ) then

          waux3(1:3) = wopos(1:3)
          waux3(2)   = 'SCALA'
          kdime      = memory_size(gevec,1_ip)
          call PAR_MAX(kdime)
          if(      waux3(3) == 'NPOIN' ) then
             dim_max = npoin
          else if( waux3(3) == 'NELEM' ) then
             dim_max = nelem
          else if( waux3(3) == 'NBOUN' ) then
             dim_max = nboun
          end if

          nullify(ggggg)
          if( INOTMASTER ) call memory_alloca(memor_dom,'GGGGG','outvar',ggggg,dim_max)

          do idime = 1,kdime
             wtype = intost(idime)
             if( idime < 10 ) then
                waux3(1) = wopos(1)(1:3)//'0'//trim(wtype(1:1))
             else if( idime < 100 ) then
                waux3(1) = wopos(1)(1:3)//trim(wtype(1:2))
             else if( idime < 1000 ) then               
                waux3(1) = wopos(1)(1:2)//trim(wtype(1:3))
             else              
                waux3(1) = trim(wtype(1:5))           
             end if
             if( IMASTER ) then
                call postpr(dusca,waux3,jttim,dutim)   
             else           
                do ipoin = 1,dim_max
                   ggggg(ipoin) = gevec(idime,ipoin)
                end do
                call postpr(ggggg,waux3,jttim,dutim) 
             end if
          end do
          if( INOTMASTER ) call memory_deallo(memor_dom,'GGGGG','outvar',ggggg)
       end if

    else if ( wopos(2) == 'SCALX' .and. kfl_ivari(1) > 0 ) then        

       if( INOTMASTER ) allocate( ggggg(npoin) )

       wauxi(2,1) = 'SCALA'
       wauxi(3,1) = wopos(3)           

       do ii = 1,2
          if( ii == 1 ) then
             wauxi(1,1) = wopos(1)(1:4)//'r'
             if( INOTMASTER ) then
                do ipoin = 1,npoin
                   ggggg(ipoin) = real(gescx(ipoin))
                end do
             end if
          else
             wauxi(1,1) = wopos(1)(1:4)//'i'           
             if( INOTMASTER ) then
                do ipoin = 1,npoin
                   ggggg(ipoin) = aimag(gescx(ipoin))
                end do
             end if
          end if
          if( INOTMASTER ) then
             call postpr(ggggg,wauxi(:,1),jttim,dutim)   
          else 
             call postpr(dusca,wauxi(:,1),jttim,dutim)   
          end if
       end do

       if( INOTMASTER ) deallocate( ggggg )

    else if ( wopos(2) == 'VECTX' .and. kfl_ivari(1) > 0 ) then        

       if( INOTMASTER ) allocate( vvvvv(ndime,npoin) )

       wauxi(2,1) = 'VECTO'
       wauxi(3,1) = wopos(3)    

       do ii = 1,2
          if( ii == 1 ) then
             wauxi(1,1) = wopos(1)(1:4)//'r'
             if( INOTMASTER ) then
                do ipoin = 1,npoin
                   do idime = 1,ndime
                      vvvvv(idime,ipoin) = real(gevex(idime,ipoin))
                   end do
                end do
             end if
          else
             wauxi(1,1) = wopos(1)(1:4)//'i'           
             if( INOTMASTER ) then
                do ipoin = 1,npoin
                   do idime = 1,ndime
                      vvvvv(idime,ipoin) = aimag(gevex(idime,ipoin))
                   end do
                end do
             end if
          end if
          if( INOTMASTER ) then
             call postpr(vvvvv,wauxi(:,1),jttim,dutim)   
          else 
             call postpr(duvec,wauxi(:,1),jttim,dutim)   
          end if
       end do
       if( INOTMASTER ) deallocate( vvvvv )

    else if ( wopos(2) == 'TENSO' ) then

       nullify(ggggg)
       if( INOTMASTER .and. kfl_ivari(1) > 0 ) call memory_alloca(memor_dom,'GGGGG','outvar',ggggg,npoin)

       if ( ndime == 1 ) then

          call runend('CANNOT POSTPROCESS A TENSOR IN 1D')

       else if ( ndime == 2 ) then        

          if( kfl_ivari(1) > 0 ) then

             waux2(1) = 'XX'        
             waux2(2) = 'YY'
             waux2(3) = 'XY'

             waux2(4) = 'TH'        
             waux2(5) = 'NU'
             waux2(6) = 'ER'

             do itens = 1,ntens
                wauxi(1,itens) = wopos(1)(1:3)//waux2(itens)
                wauxi(2,itens) = 'SCALA'
                wauxi(3,itens) = wopos(3)           
                if( INOTMASTER ) then
                   do ipoin = 1,npoin
                      ggggg(ipoin) = gevec(itens,ipoin)
                   end do
                   call postpr(ggggg,wauxi(:,itens),jttim,dutim)   
                else 
                   call postpr(dusca,wauxi(:,itens),jttim,dutim)   
                end if
             end do

          end if

       else if ( ndime == 3 ) then 

          if( kfl_ivari(1) > 0 ) then
             waux2(1) = 'XX'
             waux2(2) = 'YY'
             waux2(3) = 'ZZ'
             waux2(4) = 'YZ'
             waux2(5) = 'XZ'
             waux2(6) = 'XY'
             do itens = 1,ntens
                wauxi(1,itens) = wopos(1)(1:3)//waux2(itens)
                wauxi(2,itens) = 'SCALA'
                wauxi(3,itens) = wopos(3)           
                if( INOTMASTER ) then
                   do ipoin = 1,npoin
                      ggggg(ipoin) = gevec(itens,ipoin) 
                   end do
                   call postpr(ggggg,wauxi(:,itens),jttim,dutim)   
                else
                   call postpr(dusca,wauxi(:,itens),jttim,dutim)   
                end if
             end do
          end if
       end if

       if( INOTMASTER .and. kfl_ivari(1) > 0 ) call memory_deallo(memor_dom,'GGGGG','outvar',ggggg)

    end if

10  continue

    if( iasca == 1 ) call memgen(2_ip,one,zero)
    if( iavec == 1 ) call memgen(2_ip,one,one)
    nullify(gesca)
    nullify(gevec)
    if( if_interpolated ) then
       call interp_fe_deallocate(gevec_loc)
       call interp_fe_deallocate(gesca_loc)
    end if
    ivapo = 0

  end subroutine outvar

end module mod_outvar
!> @}
