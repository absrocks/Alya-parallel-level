subroutine nsi_elmdir(&
     order,itask,pnode,pevat,ndofn,lnods,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmdir
  ! NAME 
  !    nsi_elmdir
  ! DESCRIPTION
  !   This routine modifies the element stiffness matrix of the Navier-
  !   Stokes equations to impose the correct boundary conditions for
  !   local boundary conditions (in skew-systems). In this case the rotation 
  !   of the nodal matrices if this system is the tangent is understood to be  
  !
  !    { NORMAL , TANGENT 1 , TANGENT 2 }
  !
  !   The tangent vectors have been computed according to a good conditioning
  !   of the calculations. Therefore, only the possibility of prescribing
  !   the normal component must be considered when this option is used.
  !   Also, if the flow is confined a pressure value is prescribed.
  !   According to the value of ORDER and ITASK, the following is performed:
  !   ITASK = 1 ... Momentum + continuity equation
  !           2 ... Momentum equation
  !           3 ... Continuity equation
  !   ORDER = 1 ... Velocity and pressure
  !         = 2 ... Saving constant matrix
  !         = 3 ... This is a solid material. Put velocity and pressure 
  !                 to zero on solid nodes
  ! USES
  !    nsi_rotmat
  ! USED BY
  !    nsi_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,npoin,exnor,skcos,lpoty,nperi,&
       &                        coord,lmate,lmatn
  use def_nastin, only       :  kfl_confi_nsi,nodpr_nsi,kfl_local_nsi,&
       &                        kfl_fixno_nsi,kfl_fixpr_nsi,&
       &                        kfl_fixrs_nsi,bvess_nsi,&
       &                        skcos_nsi,kfl_algor_nsi,kfl_regim_nsi,&
       &                        kfl_perip_nsi,lperp_nsi,valpr_nsi
  implicit none
  integer(ip), intent(in)    :: order,itask,pnode,pevat,ndofn
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elmat(pevat,pevat),elrhs(pevat)
  real(rp)                   :: adiag,worma(ndime,ndime)
  integer(ip)                :: inode,ipoin,idofp,ievav,idime
  integer(ip)                :: jevav,iroty,ibopo,jnode,iperi,ipres

  if(order/=3) then
     !
     ! Prescribe one pressure degree of freedom if the flow is confined
     !
     if(kfl_confi_nsi==1.and.itask/=2) then
        do inode=1,pnode
           ipoin=lnods(inode)
           ipres=0
           if(ipoin==nodpr_nsi) then
              ipres=1
           end if
           if(ipres==1) then
              idofp=inode*ndofn
              adiag=elmat(idofp,idofp)
              do jevav=1,pevat
                 elrhs(jevav)       = elrhs(jevav)-elmat(jevav,idofp)*valpr_nsi
                 elmat(jevav,idofp) = 0.0_rp 
                 elmat(idofp,jevav) = 0.0_rp
              end do
              if(order==1) then
                 elmat(idofp,idofp) = adiag
                 elrhs(idofp)       = valpr_nsi*adiag
              end if
           end if
        end do
        if(nperi>0) then
           do inode=1,pnode
              ipoin=lnods(inode)
              ipres=0
              if(kfl_perip_nsi==1) then
                 iperi=1
                 loop_iperi: do while(lperp_nsi(iperi)/=0)
                    if(lperp_nsi(iperi)==ipoin) then
                       ipres=1
                       exit loop_iperi
                    end if
                    iperi=iperi+1
                 end do loop_iperi
              else
                 if(ipoin==nodpr_nsi) ipres=1
              end if
              if(ipres==1) then
                 idofp=inode*ndofn
                 adiag=elmat(idofp,idofp)
                 do jevav=1,pevat
                    elrhs(jevav)       = elrhs(jevav)-elmat(jevav,idofp)*valpr_nsi
                    elmat(jevav,idofp) = 0.0_rp 
                    elmat(idofp,jevav) = 0.0_rp
                 end do
                 if(order==1) then
                    elmat(idofp,idofp) = adiag
                    elrhs(idofp)       = valpr_nsi*adiag
                 end if
              end if
           end do
        else
           do inode=1,pnode
              ipoin=lnods(inode)
              if(ipoin==nodpr_nsi) then
                 idofp=inode*ndofn
                 adiag=elmat(idofp,idofp)
                 if(adiag==0.0_rp) adiag=1.0_rp
                 do jevav=1,pevat
                    elrhs(jevav)       = elrhs(jevav)-elmat(jevav,idofp)*valpr_nsi
                    elmat(jevav,idofp) = 0.0_rp
                    elmat(idofp,jevav) = 0.0_rp
                 end do
                 if(order==1) then
                    elmat(idofp,idofp) = adiag
                    elrhs(idofp)       = valpr_nsi*adiag
                 end if
              end if
           end do
        end if
     end if
     !
     ! Rotate the nodal matrices, the element nodal velocities and RHS to 
     ! prescribe boundary conditions in a skew-system, either the tangent
     ! one or another prescribed by the user. Also, periodical boundary 
     ! conditions are accounted for.
     !
     if(itask/=3.and.itask/=5) then
        if(kfl_local_nsi==1) then
           do inode=1,pnode
              ipoin=lnods(inode)
              ibopo=lpoty(ipoin)
              if(ibopo>0) then
                 iroty=kfl_fixrs_nsi(ibopo)
                 if(iroty==-1) then                                    ! Tangent system
                    call nsi_rotmat(&
                         &      itask,inode,pnode,ndofn,pevat,&
                         &      elmat,elrhs,exnor(1,1,ibopo),worma)
                 else if(iroty>=1) then                                ! Given system
                    call nsi_rotmat(&
                         &      itask,inode,pnode,ndofn,pevat,&
                         &      elmat,elrhs,skcos(1,1,iroty),worma)
                 else if(iroty==-2) then                               ! Given system
                    call nsi_rotmat(&
                         &      itask,inode,pnode,ndofn,pevat,&
                         &      elmat,elrhs,skcos_nsi(1,1,ibopo),worma)
                 else if(iroty==-3) then                               ! Geometrical normals
                    call nsi_rotmat(&
                         &      itask,inode,pnode,ndofn,pevat,&
                         &      elmat,elrhs,skcos(1,1,ibopo),worma)
                 end if
              end if
           end do
        end if
        !  
        ! Velocity Dirichlet boundary conditions
        !
        do inode=1,pnode
           ievav=(inode-1)*ndofn
           ipoin=lnods(inode)
           do idime=1,ndime
              ievav=ievav+1
              if( kfl_fixno_nsi(idime,ipoin) > 0 ) then
                 adiag=elmat(ievav,ievav)
                 do jevav=1,pevat
                    elmat(ievav,jevav)=0.0_rp
                 end do
                 if(order==1) then 
                    elmat(ievav,ievav)=adiag
                    elrhs(ievav)=adiag*bvess_nsi(idime,ipoin,1)
                 end if
              end if
           end do
        end do
     end if

  else
     !
     ! Solid element
     !
     do ievav=1,pevat
        do jevav=1,pevat
           elmat(ievav,jevav)=0.0_rp
           elmat(jevav,ievav)=0.0_rp
        end do
        elrhs(ievav)=0.0_rp
     end do
     do inode=1,pnode
        ievav=(inode-1)*ndofn
        ipoin=lnods(inode)
        !if(lmatn(ipoin)==-1) then
        !   do idime=1,ndofn
        !      ievav=ievav+1
        !      elmat(ievav,ievav)=1.0_rp
        !   end do
        !end if
     end do

  end if

end subroutine nsi_elmdir
