!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_volume_bound.f90
!> @author  Mariano Vazquez 
!> @date    16/NOV/1966
!> @brief   Compute cavities volume
!> @details Compute cavities volume
!>           - itask=1 The end of an internal iteration
!>           - itask=2 The end of an external iteration
!> @}
!------------------------------------------------------------------------
subroutine sld_volume_bound(itask)
  use def_kintyp, only       :  ip,rp
  use def_domain
  use def_solidz
  use def_master
  implicit none

  integer(ip), intent(in) :: itask             !< 1, inner iterations; 2, end of inner iterations

  integer(ip) :: ivolu,icavi,ielem,ipoin,idime,iboun,igaub,inodb,idofn,jdime,kdime
  integer(ip) :: pelty,pblty,pnodb,pgaub,pnode,inode,ipass

  real(rp) :: xcent,ycent,zcent,node1(3),node2(3),node3(3),ax,ay,az,bx,by,bz
  real(rp) :: gpele(3),norme,venor(3),super,volum(mcavi_sld),fonct,bidon
  real(rp) :: dummr(4),super_total,deltv,deltp,compl
  real(rp) :: coord_vol(3,3)
  real(rp) :: baloc(ndime,ndime),gppush(ndime,ndime)
  real(rp) :: bocod(ndime,mnodb),elcod(ndime,mnode)
  real(rp)    :: gbsur,eucta,tract(ndime),gppre
  real(rp) :: elpush(ndime,ndime,mnodb),unr,xnr,fnr

  !
  ! https://www.sangakoo.com/en/unit/volumes-calculation-using-gauss-theorem
  !


  volum = 0.0_rp

  if( INOTMASTER ) then

     super_total = 0.0_rp
     super       = 0.0_rp

     do iboun=1,nboun
        ivolu= 0
        do icavi= 1,mcavi_sld
           if (lbset(iboun) == iocav_sld(icavi)) then
              ivolu= icavi
           end if
        end do

        if ((ivolu > 0)) then
           !
           ! Element properties and dimensions
           !
           pblty=ltypb(iboun) 
           pnodb=nnode(pblty)
           pgaub=ngaus(pblty)
           ielem=lelbo(iboun)
           pelty=ltype(ielem)
           pnode=nnode(pelty)

           !
           ! Gather operations: BOCOD, ELCOD, ELPUSH
           !
           do inodb=1,pnodb
              ipoin=lnodb(inodb,iboun)
              do idime=1,ndime
                 bocod(idime,inodb) = coord(idime,ipoin) + displ(idime,ipoin,1)
              end do
           end do
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin) + displ(idime,ipoin,1)
              end do
           end do
           !
           ! Loop over Gauss points
           !
           gauss_points: do igaub=1,pgaub
              call bouder(pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),bocod,baloc,eucta)    
              gbsur=elmar(pblty)%weigp(igaub)*eucta 
              call chenor(pnode,baloc,bocod,elcod)        ! Check normal

              gpele= 0.0_rp
              do inodb = 1,pnodb
                 do idime= 1,ndime
                    gpele(idime) = gpele(idime) + elmar(pblty)%shape(inodb,igaub) * &
                     (bocod(idime,inodb) - ocavi_sld(idime,ivolu))
                 end do
              end do

              do idime = 1,ndime
                 ! negative because normal is exterior
                 volum(ivolu) = volum(ivolu) - gpele(idime) * baloc(idime,ndime)* gbsur / 3.0_rp
              end do              

           end do gauss_points

        end if

     end do !iboun

  end if !NOTMASTER

  call pararr('SUM',0_ip,mcavi_sld,volum)

  call pararr('SUM',0_ip,1_ip,super_total)

  do ivolu=1,mcavi_sld !loop over the cavities, cavities defined in solidz module

     volst(ITER_K,ivolu) =   volum(ivolu)

     !update volst, only at the end of the iterations loop
     if (itask == 2) then
        volst(TIME_N_MINUS_2,ivolu) = volst(TIME_N_MINUS_1,ivolu)  !V^n-2
        volst(TIME_N_MINUS_1,ivolu) = volst(TIME_N,ivolu)          !V^n-1
        volst(TIME_N,ivolu)         = volst(ITER_K,ivolu)          !V^n
     end if
     if (ittim == 1) then
        !volvr_sld(1) = volst(TIME_N,1)
        !volvr_sld(2) = volst(TIME_N,2)
        volst(TIME_N,ivolu)         = volst(ITER_K,ivolu)          !V^n
        volst(TIME_N_MINUS_2,ivolu) = volst(TIME_N,ivolu)
        volst(TIME_N_MINUS_1,ivolu) = volst(TIME_N,ivolu)  !initially V^n-1 = V^n = V^i 
     end if
     
     if (mcavi_sld == 1) then
        !
        ! When there is one ventricle only, put 1.0 in the second one volume
        !
        volst(ITER_K,2)=1.0_rp

        !deltaV
        volst(TIME_N,2) = 1.0_rp

     end if

     !
     ! Broadcast values to the slaves: everyone need to know the current volume value for windkessel models
     !  
     call pararr('BCT',0_ip,12_ip,volst(ITER_K,ivolu))
     call pararr('BCT',0_ip,12_ip,volst(TIME_N,ivolu))

  end do

100 format (2(F16.8,','))  
101 format (6(F19.11,',')) 
end subroutine sld_volume_bound
