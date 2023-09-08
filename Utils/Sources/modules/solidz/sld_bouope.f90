!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_bouope.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine computes boundary operations
!> @details This routine computes boundary operations
!> @}
subroutine sld_bouope()
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_bouope
  ! NAME
  !    sld_bouope
  ! DESCRIPTION
  !    boundary operations:
  !    1. Compute elemental matrix and RHS
  !    2. Impose Dirichlet boundary conditions
  !    3. Assemble them
  ! USES
  !    sld_bouave
  !    sld_bougat
  !    sld_elmpro
  !    bouder
  !    chenor
  !    sld_bouwal
  !    elmder
  !    cartbo
  !    sld_bouopb
  !    sld_elmdir
  !    sld_assrhs
  !    assmat
  ! USED BY
  !    sld_matrix
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_solidz
  !use mod_sld_commdom, only: commdom_sld_space_funcion

  implicit none

  integer(ip) :: ielem,ipoin,idime,jdime,iboun,igaub,inodb,idofn,icavi
  integer(ip) :: pelty,pblty,pnodb,pgaub,pnode,inode,kfixno,itott
  real(rp)    :: elrhs(ndofn_sld*mnode),elfex(ndofn_sld*mnode)
  real(rp)    :: baloc(ndime,ndime),gppush(ndime,ndime)
  real(rp)    :: bocod(ndime,mnodb),elcod(ndime,mnode),elpre(mnodb),elpush(ndime,ndime,mnodb),unr,xnr,fnr
  real(rp)    :: gbsur,eucta,tract(ndime),tract_ref(ndime),kspring,muvisco,gppre,unormal,xnormal,fnormal,tnormal

  !
  ! Loop over boundaries
  !
  !  call commdom_sld_space_funcion( bou_prop=bvnat_sld(1,1:nboun,1), fixbo=7_ip, idime=1_ip )  !< 2015Jul30

  if (kfl_icodb .ne. 0) then

     boundaries: do iboun=1,nboun

        if(kfl_fixbo_sld(iboun)==2 &
             .or. kfl_fixbo_sld(iboun)==3 &
             .or. kfl_fixbo_sld(iboun)==6 &
             .or. kfl_fixbo_sld(iboun)==7 &                                                  !< 2015Jul30
             .or. kfl_fixbo_sld(iboun)==8 &
             .or. kfl_fixbo_sld(iboun)==9 &
             .or. ( kfl_fixbo_sld(iboun)==4 .and. kfl_cycle_sld==1 ) ) then
           !
           ! fixbo=2  --> Nodal pressure (in press)
           ! fixbo=3  --> Boundary element pressure (in bvnat, as a traction on a LOCAL basis)
           ! fixbo=5  --> Nodal traction coming from a field (in bvess)
           ! fixbo=6  --> Boundary element traction (in bvnat, as a traction on a GLOBAL basis)
           ! fixbo=8  --> Damper-like boundary condition
           ! fixbo=9  --> Boundary element Windkessel pressure (in bvnat, as a traction on a LOCAL basis)
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
           ! Initialize
           !
           elrhs(:) = 0.0_rp
           elfex(:) = 0.0_rp
           !
           ! Gather operations: BOCOD, ELCOD, ELPUSH
           !

           do inodb=1,pnodb
              ipoin=lnodb(inodb,iboun)
              do idime=1,ndime
                 bocod(idime,inodb) = coord(idime,ipoin)
                 elpre(inodb)       = press(ipoin,1_ip)
                 if( kfl_gdepo /= 0 ) then
                    do jdime=1,ndime
                       elpush(idime,jdime,inodb)= gdeinv(idime,jdime,ipoin)*gdedet(ipoin)
                    end do
                 end if
              end do
           end do
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              do idime=1,ndime
                 elcod(idime,inode)=coord(idime,ipoin)! + displ(idime,ipoin,1)
              end do
           end do
           !
           ! Loop over Gauss points
           !
           gauss_points: do igaub=1,pgaub
              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&
                   bocod,baloc,eucta)
              gbsur=elmar(pblty)%weigp(igaub)*eucta
              call chenor(&                                           ! Check normal
                   pnode,baloc,bocod,elcod)

              !
              ! Compute push forward operator in the gauss point: F^{-T}*J
              !
              gppush= 0.0_rp
              gppre= 0.0_rp
              if( kfl_gdepo /= 0 ) then
                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    gppre = gppre + elmar(pblty)%shape(inodb,igaub) * elpre(inodb)
                    do idime=1,ndime
                       do jdime=1,ndime
                          gppush(idime,jdime) = gppush(idime,jdime) + elmar(pblty)%shape(inodb,igaub) * elpush(idime,jdime,inodb)
                       end do
                    end do
                 end do
              end if

              if(kfl_fixbo_sld(iboun)==2) then
                 !
                 ! pressure force(igaus) = n_i * p(igaus) * gbsur
                 !                       = baloc(1:ndime,ndime) * p(igaus) * gbsur
                 !
                 ! baloc(1:ndime,ndime) --> exterior normal vector
                 !
                 gppre= 0.0_rp
                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    gppre = gppre + elmar(pblty)%shape(inodb,igaub) * elpre(inodb)
                 end do

                 tract(1:ndime) = baloc(1:ndime,ndime) * gppre

              else if (kfl_fixbo_sld(iboun)==3) then

                 ! 1_ip  means "boundary elements"
                 if (nfunc_sld > 0) call sld_funbou( 1_ip, iboun,bvnat_sld(1,iboun,1))

                 gppre = bvnat_sld(1,iboun,1)

                 tract(1:ndime) = baloc(1:ndime,ndime) * gppre

              else if (kfl_fixbo_sld(iboun)==4 .and. kfl_cycle_sld==1) then
                 !The CCMT is called from sld_endite
                 if (kfl_codbo(iboun)==1) then
                    gppre=ptota_sld(1)  !LV
                 else
                    gppre=ptota_sld(1)/2.0_rp  !RV: simply scale the pressure calculated bu the LV
                    ptota_sld(2)=gppre !for writing
                    !gppre=0.0_rp
                 end if

                 tract(1:ndime) = baloc(1:ndime,ndime) * gppre

              else if (kfl_fixbo_sld(iboun)==6) then

                 tract(1:ndime) = bvnat_sld(1:ndime,iboun,1)

              else if (kfl_fixbo_sld(iboun)==7) then                               !< 2015Jul30
                 !
                 !  bvnat_sld(1,iboun,1) == F(x)
                 !
                 tract(1:ndime) = baloc(1:ndime,ndime) * bvnat_sld(1,iboun,1)

              else if (kfl_fixbo_sld(iboun)==8) then
                 !
                 !  normal damper-like forces
                 !  F = -k x + mu dx/dt
                 unormal= 0.0_rp
                 xnormal= 0.0_rp
                 fnormal= 0.0_rp
                 kfixno = 0_ip
                 do inodb = 1,pnodb
                    unr= 0.0_rp
                    xnr= 0.0_rp
                    fnr= 0.0_rp
                    ipoin = lnodb(inodb,iboun)
                    itott = (ipoin             -1)*ndime
                    idofn=  (lboel(inodb,iboun)-1)*ndime
                    if (kfl_fixno_sld(1,ipoin) == 2) kfixno = 2
                    do idime=1,ndime
                       idofn=idofn+1
                       itott=itott+1
                       fnr= fnr + baloc(idime,ndime)*rhsid(itott) ! forces nodal normal projection
                       unr= unr + baloc(idime,ndime)*veloc_sld(idime,ipoin,ITER_K)  ! velocity nodal normal projection
                       xnr= xnr + baloc(idime,ndime)*displ(idime,ipoin,ITER_K)  ! displacement nodal normal projection
                    end do
                    unormal = unormal + elmar(pblty)%shape(inodb,igaub) * unr ! velocity normal projection at the gauss point
                    xnormal = xnormal + elmar(pblty)%shape(inodb,igaub) * xnr ! displacement normal projection at the gauss point
                    fnormal = fnormal + elmar(pblty)%shape(inodb,igaub) * fnr ! force normal projection at the gauss point
                 end do

                 ! tests:
                 ! if (fixed normal displacement for at least one node) then check fnormal
                 !    if (fnormal < 0.0) then
                 !       free node normal displacement
                 !       compute tnormal and add it to elrhs
                 !    else
                 !       keep node normal displacement fixed
                 !       keep tnormal = 0
                 !    end if
                 ! else
                 !    if (xnormal < 0.0) then
                 !       keep node normal displacement free
                 !       compute tnormal and add it to elrhs
                 !    else
                 !       fix node normal displacement
                 !       keep tnormal = 0
                 !    end if
                 ! end if

                 kspring= - bvnat_sld(1,iboun,1)
                 muvisco= - bvnat_sld(2,iboun,1)
                 tnormal= 0.0_rp

                 if (kfixno == 2) then
                    !if (fnormal .le. 0.0_rp) then
                       kfixno = 0
                       tnormal = kspring*xnormal + muvisco * abs(unormal) !normal traction component
                    !end if
                 else
                    !if (xnormal .le. 0.0_rp) then
                       kfixno = 0
                       tnormal = kspring*xnormal + muvisco * abs(unormal) !normal traction component
                    !else
                    !   kfixno = 2
                    !end if
                 end if

                 do inodb=1,pnodb
                    ipoin=lnodb(inodb,iboun)
                    if(kfl_fixno_sld(1,ipoin)/=1) then
                       kfl_fixno_sld(1,ipoin) = kfixno
                    end if
                 end do

                 ! tnormal = 100000.0_rp   ! debug
                 !! leave this commented: it is to check that all is ok: tnormal positive means force outwards (i.e. following normal)

                 !
                 ! if the node moves against the normal (i.e. inwards), xnormal is negative. then
                 ! a kspring negative produces a tnormal positive.
                 ! then, a force in the normal direction is generated to compensate movement
                 !
                 ! total traction, recall that baloc(1:ndime,ndime) is the normal versor
                 tract(1:ndime) = baloc(1:ndime,ndime) * tnormal

              else if (kfl_fixbo_sld(iboun)==9) then

                 icavi= int(bvnat_sld(1,iboun,1))
                 tract(1:ndime) = - baloc(1:ndime,ndime) * pwink(ITER_K,icavi)     ! minus, because normal is defined pointing outwards

              end if  !kfl_fixbo_sld(iboun)

              !
              ! Push forward tract vector: tract= tract_ref * J * F^{-T}
              !

 ! This was a bug. Push Forward is actually not required.
 !
 !             if (kfl_cshel_sld(1) == 0 .and. kfl_cshel_sld(2) == 0 .and. kfl_cshel_sld(3) == 0) then        ! Only do this when there are no continuum shell elements
 !
 !                if ((kfl_fixbo_sld(iboun) == 3) .or. (kfl_fixbo_sld(iboun) == 9)) then
 !                   tract_ref = tract
 !                   tract     = 0.0_rp
 !                   do idime= 1,ndime
 !                      do jdime= 1,ndime
 !                         tract(idime)= tract(idime)+gppush(jdime,idime)*tract_ref(jdime)
 !                      end do
 !                   end do
 !                end if
 !             end if

              do inodb=1,pnodb
                 idofn=(lboel(inodb,iboun)-1)*ndime
                 ipoin=lnodb(inodb,iboun)
                 do idime=1,ndime
                    idofn=idofn+1
                    if(kfl_fixno_sld(idime,ipoin)/=1) then
                       if(kfl_fixno_sld(idime,ipoin)/=2) then
                          if(kfl_fixno_sld(idime,ipoin)/=3) then
                             elrhs(idofn)=elrhs(idofn) + gbsur*tract(idime) * elmar(pblty)%shape(inodb,igaub)
                             !
                             ! External forces with Newmann part
                             elfex(idofn)=elfex(idofn) + gbsur*tract(idime) * elmar(pblty)%shape(inodb,igaub)
                          end if
                       end if
                    end if
                 end do
              end do

           end do gauss_points
           !
           ! Assembly
           !
           call assrhs(&
                ndofn_sld,pnode,lnods(1,ielem),elrhs,rhsid)
           call assrhs(&
                ndofn_sld,pnode,lnods(1,ielem),elfex,fexte_sld)

        end if  !kfl_fixbo_sld(iboun)==2 .or. 3 .or. kfl_codbo(iboun) /= 100

     end do boundaries

  end if

end subroutine sld_bouope
