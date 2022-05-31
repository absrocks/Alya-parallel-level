!------------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_elerhs_matrix.f90
!> @date    01/08/2012
!> @author  Mariano Vazquez
!> @brief   Computes the elrhs matrix and vector
!> @details Computes the elrhs matrix and vector
!> @}
!------------------------------------------------------------------------
subroutine nsa_elerhs_matrix(&
     elrhs,elsax,elmat,taudi,pnode,pevat,inode,igaus,ndaux,&
     dvolu,elunk,xvelo,xunkn,&
     xconv,xdiff,dconv,ddiff,&
     xconv_der,&
     gunkn,xshap,cartd,hesma,&
     xtide,xlopr,&
     xsube,shote,shmet,ielem)
  use      def_master
  use      def_domain
  use      def_nastal
  implicit none
  integer(ip) :: &
       pnode,pevat,ndaux,ndofi,jdofn,idofn,kdofn,idime,jdime,kdime,ievat,jevat,inode,jnode,igaus,ielem,&
       iauxi,ipoin,itott,ITER_NEWTON,if_assemble_rhs, if_assemble_mat

  real(rp)    :: &
       dvolu,elrhs(pevat),elsax(pevat),elmat(pevat,pevat),&
       state(ndofn_nsa)  ,galte(ndofn_nsa)  ,shote(ndofn_nsa),rhsmat(ndofn_nsa),&
       state_n(ndofn_nsa),galte_n(ndofn_nsa),shote_n(ndofn_nsa),&
       aumat,autim(ndofn_nsa),autim_n(ndofn_nsa),ausax(ndofn_nsa),&
       xconv(ndofn_nsa,ndofn_nsa,ndime),elunk(ndofn_nsa,mnode,ncomp_nsa),&
       xconv_der(ndofn_nsa,ndofn_nsa,ndofn_nsa,ndime),&
       xconv_newt(ndofn_nsa,ndofn_nsa,ndime),&
       gunkn(ndofn_nsa,ndime),xshap(mnode),cartd(ndime,mnode),&
       ddiff(ndofn_nsa,ndofn_nsa,2,ndime), &
       xsube(ndofn_nsa,mgaus,3),xdiff(ndofn_nsa,ndofn_nsa,ndime,ndime),&
       hesma(ndime,ndime,mnode), &
       xunkn(ndofn_nsa,mgaus,3),xvelo(ndime),shmet(ndime,ndime,ndofn_nsa),&
       dconv(ndofn_nsa,ndofn_nsa),profi(ndofn_nsa,ndofn_nsa),taudi(ndofn_nsa),diff_fact,supg_fact,&
       ximpl_visc, ximpl_conv, ximpl_stab, ximpl_shot

  real(rp) :: &
       advec_matrix(ndofn_nsa,ndofn_nsa,mnode),& ! Advective part
       diffu_matrix(ndofn_nsa,ndofn_nsa,mnode),& ! Diffusive part
       shote_matrix(ndofn_nsa,ndofn_nsa,mnode),& ! Shock capturing part
       debudiagi(ndofn_nsa,ndofn_nsa),&
       stabi_matrix(ndofn_nsa,ndofn_nsa,mnode)  ,& ! Stabilization part = Adjoint matrix * Diagonal subscale part
       timas_matrix(ndofn_nsa,ndofn_nsa,mnode,2),& ! Temporal part for time and pseudotime 
       adjoi_matrix(ndofn_nsa,ndofn_nsa),&       ! Adjoint matrix
       subdi_matrix(ndofn_nsa,ndofn_nsa,mnode),& ! Diagonal subscale part
       subes(ndofn_nsa),&                        ! Subscale
       xlopr(ndofn_nsa,ndofn_nsa),&              ! 
       xtide(ndofn_nsa),&              ! 
       auele_rhs,auele(pevat),elunk_jevat,resid_n, dtinv_eqs(5,2),elunk_value

  real(rp) :: xx_taudi,cn_left,cn_pseudo,explicit_pseudo(5)
    
  cn_left= 1.0_rp  
  if (kfl_tisch_nsa == 3) cn_left= 0.5_rp  ! crank-nicolson 0.5
  cn_pseudo= 1.0_rp - cn_left

  ximpl_visc= 0.0_rp
  if (kfl_ximpl_nsa(1)==1) ximpl_visc= 0.5_rp
  ximpl_conv= 0.0_rp
  if (kfl_ximpl_nsa(2)==1) ximpl_conv= 0.5_rp
  ximpl_stab= 0.0_rp
  if (kfl_ximpl_nsa(3)==1) ximpl_stab= 0.5_rp
  ximpl_shot= 0.0_rp
  if (kfl_ximpl_nsa(4)==1) ximpl_shot= 0.5_rp
  
  supg_fact= 1.0_rp
  if (kfl_stabi_nsa == 6) supg_fact=0.0_rp  !generalized supg

  ITER_NEWTON= ITER_K

  ipoin= lnods(inode,ielem)

  !
  ! PHYSICAl time step:
  ! kfl_dttyp_nsa defines if it is local or not 
  ! it is initialized with dtinv_nsa
  !
  dtinv_eqs(1:5,DT_PHYSICAL)    = dtinv_nsa
  if (kfl_dttyp_nsa(1) > 0 ) then    ! momentum, local time step
     do idime=1,ndime
        itott= (ipoin-1) * ndofn_nsa + idime
        dtinv_eqs(idime,DT_PHYSICAL)= 1.0_rp/dtieq_nsa(1,ipoin,DT_PHYSICAL) 
     end do
  end if  
  if (kfl_dttyp_nsa(2) > 0 ) then    ! continuity, local time step
     itott= (ipoin-1) * ndofn_nsa + ndime + 1
     dtinv_eqs(ndime+1,DT_PHYSICAL)= 1.0_rp/dtieq_nsa(2,ipoin,DT_PHYSICAL) 
  end if  
  if (kfl_dttyp_nsa(3) > 0 ) then    ! energy, local time step
     itott= (ipoin-1) * ndofn_nsa + ndime + 2
     dtinv_eqs(ndime+2,DT_PHYSICAL)= 1.0_rp/dtieq_nsa(3,ipoin,DT_PHYSICAL) 
  end if  

  !
  ! PSEUDO time step:
  ! being non-physical, it is always local
  !
!  dtinv_eqs(1:ndime,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(1,ipoin,DT_PHYSICAL) * safet_nsa / safet_pseud_nsa
!  dtinv_eqs(ndime+1,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(2,ipoin,DT_PHYSICAL) * safet_nsa / safet_pseud_nsa
!  dtinv_eqs(ndime+2,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(3,ipoin,DT_PHYSICAL) * safet_nsa / safet_pseud_nsa

  if (kfl_pseud_nsa == 1) then
     dtinv_eqs(1:ndime,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(1,ipoin,DT_PSEUDO)
     dtinv_eqs(ndime+1,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(2,ipoin,DT_PSEUDO)
     dtinv_eqs(ndime+2,DT_PSEUDO)  = 1.0_rp/dtieq_nsa(3,ipoin,DT_PSEUDO)
  else
     dtinv_eqs(1:ndime,DT_PSEUDO)  = 0.0_rp
     dtinv_eqs(ndime+1,DT_PSEUDO)  = 0.0_rp
     dtinv_eqs(ndime+2,DT_PSEUDO)  = 0.0_rp
  end if

  explicit_pseudo= 1.0_rp
  if (kfl_pseud_nsa == 1) then
     ! correct to use the pseudo-time safety factor, which can be different than the time safety factor
!!!     dtinv_eqs(1:5,DT_PSEUDO) = dtinv_eqs(1:5,DT_PHYSICAL) * safet_nsa / safet_pseud_nsa
!     if (itinn(modul) == 1) then
!        dtinv_eqs(1:5,DT_PSEUDO) = dtinv_eqs(1:5,DT_PHYSICAL)
!        dtinv_eqs(1:5,DT_PHYSICAL) = 0.0_rp
!     end if
  end if
  

  do idofn=1,ndofn_nsa
     if (itinn(modul) > 1) then
        if (kfl_pseud_nsa == 1) then        
           explicit_pseudo(idofn) = dtinv_eqs(idofn,DT_PSEUDO) + dtinv_eqs(idofn,DT_PHYSICAL) 
           if (explicit_pseudo(idofn) > 0.0_rp) then
              explicit_pseudo(idofn)= dtinv_eqs(idofn,DT_PHYSICAL) / explicit_pseudo(idofn)
           end if
        end if
     end if

     state(idofn)= 0.0_rp
     shote(idofn)= 0.0_rp    ! shote is recomputed from shmet matrix
     galte(idofn)= 0.0_rp
     autim(idofn)= 0.0_rp     
     rhsmat(idofn)= 0.0_rp     
     state_n(idofn)= 0.0_rp
     galte_n(idofn)= 0.0_rp
     shote_n(idofn)= 0.0_rp
     autim_n(idofn)= 0.0_rp     
     subes(idofn)= 0.0_rp
     ausax(idofn)= 0.0_rp     
     do jdofn=1,ndaux  
        adjoi_matrix(jdofn,idofn) = 0.0_rp ! Adjoint matrix
        xconv_newt(jdofn,idofn,1:ndime) = 0.0_rp 
        do jnode=1,pnode
           advec_matrix(idofn,jdofn,jnode) = 0.0_rp ! Advective part
           diffu_matrix(idofn,jdofn,jnode) = 0.0_rp ! Diffusive part
           shote_matrix(idofn,jdofn,jnode) = 0.0_rp ! Shock capturing part
           subdi_matrix(idofn,jdofn,jnode) = 0.0_rp ! Diagonal subscale part
           stabi_matrix(idofn,jdofn,jnode) = 0.0_rp ! Stabi part = Adjoint matrix * Diagonal subscale part
           timas_matrix(idofn,jdofn,jnode,DT_PSEUDO) = 0.0_rp ! Temporal part
           timas_matrix(idofn,jdofn,jnode,DT_PHYSICAL)   = 0.0_rp ! Temporal part
        end do
     enddo
  enddo

  if (kfl_linea_nsa == 2) then
     !
     ! Compute NR contribution
     !
     do idofn= 1,ndaux
        do jdofn= 1,ndaux
           do kdofn= 1, ndaux              
              xconv_newt(idofn,jdofn,1:ndime)= xconv_newt(idofn,jdofn,1:ndime) &
                   + xconv_der(idofn,jdofn,kdofn,1:ndime) * (xunkn(kdofn,igaus,ITER_NEWTON) - xunkn(kdofn,igaus,TIME_N))
           end do
        end do
     end do     
  end if

  do idime=1,ndime
     do idofn=1,ndofn_nsa
        do jdofn=1,ndaux  
           xx_taudi = taudi(jdofn)
           do jnode=1,pnode
              ! Advective part
              advec_matrix(idofn,jdofn,jnode) = advec_matrix(idofn,jdofn,jnode) &
                   + (xconv(idofn,jdofn,idime) + xconv_newt(idofn,jdofn,idime)) &
                   * cartd(idime,jnode) * xshap(inode)
              ! Diagonal subscale part
              diff_fact= supg_fact * ddiff(idofn,jdofn,1,idime)
              subdi_matrix(idofn,jdofn,jnode) = subdi_matrix(idofn,jdofn,jnode) &
                   + xconv(idofn,jdofn,idime) * cartd(idime,jnode) * xx_taudi &
                   - diff_fact * cartd(idime,jnode) * xx_taudi

              do jdime = 1,ndime
                 diff_fact= supg_fact *xdiff(idofn,jdofn,idime,jdime) 
                 subdi_matrix(idofn,jdofn,jnode) = subdi_matrix(idofn,jdofn,jnode) &
                      - diff_fact * hesma(idime,jdime,jnode) * xx_taudi
                 ! Diffusive part
                 diffu_matrix(idofn,jdofn,jnode) = diffu_matrix(idofn,jdofn,jnode) &
                      + xdiff(idofn,jdofn,idime,jdime) * cartd(jdime,jnode) * cartd(idime,inode)
              end do

              if (idofn == jdofn) then
                 do jdime = 1,ndime
                    ! Shock capturing
                    shote_matrix(idofn,jdofn,jnode) = shote_matrix(idofn,jdofn,jnode) &
                         + shmet(jdime,idime,idofn) * cartd(jdime,jnode) * cartd(idime,inode)
                 end do
              end if

           end do
        end do        
     end do
  end do

!!!!! ojooooooo, es para hacer pruebis
!!!!  advec_matrix= 1.0_rp


  do idofn=1,ndofn_nsa        
     do jdofn=1,ndaux          
        ! Adjoint matrix        
        if (kfl_lopre_nsa == 0) then ! (MM) This 'if' should disappear in the future
           adjoi_matrix(jdofn,idofn) = dconv(jdofn,idofn) * xshap(inode)
        else
           adjoi_matrix(jdofn,idofn) = 0.0_rp           
        end if
        do idime=1,ndime
           diff_fact= supg_fact *ddiff(jdofn,idofn,2,idime)
           adjoi_matrix(jdofn,idofn) = adjoi_matrix(jdofn,idofn) &
                + (xconv(jdofn,idofn,idime) + diff_fact) * cartd(idime,inode)

           do jdime = 1,ndime
              diff_fact= supg_fact *xdiff(idofn,jdofn,idime,jdime)
              adjoi_matrix(idofn,jdofn) = adjoi_matrix(idofn,jdofn) &
                  + diff_fact * hesma(idime,jdime,inode)
           end do
           
        end do        
     end do
  end do

!!$
  do jnode=1,pnode
     do idofn= 1,ndofn_nsa
        ! Temporal part
        timas_matrix(idofn,idofn,jnode,DT_PSEUDO)    = xshap(inode) * xshap(jnode) * dtinv_eqs(idofn,DT_PSEUDO)

        if (kfl_lopre_nsa > 0) then  
           !
           ! If preconditioning and pseudo-time, then transform physical time mass matrix. 
           !  System matrix comes transformed from nsa_gauvalxy.
           !
           do jdofn= 1,ndofn_nsa
              timas_matrix(idofn,jdofn,jnode,DT_PHYSICAL)  = &
                   xshap(inode) * xshap(jnode) * dtinv_eqs(idofn,DT_PHYSICAL) * xlopr(idofn,jdofn)
           end do
        else
           timas_matrix(idofn,idofn,jnode,DT_PHYSICAL)  = xshap(inode) * xshap(jnode) * dtinv_eqs(idofn,DT_PHYSICAL)
        end if
        
        if (kfl_stabi_nsa >= 1) then
           do jdofn= 1,ndofn_nsa
              do kdofn= 1,ndofn_nsa
                 ! Stabilization part = Adjoint matrix * Diagonal subscale part
                 stabi_matrix(idofn,jdofn,jnode)= stabi_matrix(idofn,jdofn,jnode) &
                      + adjoi_matrix(idofn,kdofn) * subdi_matrix(kdofn,jdofn,jnode)
              end do
           end do
        end if
     end do
  end do
  
  do idofn= 1,ndaux
     ievat = (inode-1) * ndofn_nsa + idofn
     elsax(ievat) = 0.0_rp
     
     do jdofn=1,ndaux
        do jnode=1,pnode
           if_assemble_rhs= 1

           if (kfl_linea_nsa == 3) then
              if_assemble_rhs = 0
              if (itinn(modul) > 1) then
                 if (jnode == inode) then
                    if (jdofn == idofn) if_assemble_rhs = 1
                 end if
              end if
           end if
           jevat= (jnode-1) * ndofn_nsa + jdofn

           ! explicit formulation terms
           galte(idofn)= galte(idofn) &
                + (advec_matrix(idofn,jdofn,jnode) &
                +  diffu_matrix(idofn,jdofn,jnode)) * elunk(jdofn,jnode,ITER_NEWTON)

           state(idofn)= state(idofn) &
                + (stabi_matrix(idofn,jdofn,jnode))* elunk(jdofn,jnode,ITER_NEWTON) 

           shote(idofn)= shote(idofn) &
                + (shote_matrix(idofn,jdofn,jnode))* elunk(jdofn,jnode,ITER_NEWTON) 
           subes(jdofn)= subes(jdofn) - subdi_matrix(jdofn,idofn,jnode) * elunk(idofn,jnode,ITER_NEWTON)
           
           if_assemble_mat = 1
           if (kfl_linea_nsa == 3) then
              if (if_assemble_rhs == 1) if_assemble_mat = 0
           end if
           aumat= 0.0_rp
           if (if_assemble_mat == 1) then
               aumat = &
                   advec_matrix(idofn,jdofn,jnode) &
                   + diffu_matrix(idofn,jdofn,jnode) &
                   + stabi_matrix(idofn,jdofn,jnode) &
                   + shote_matrix(idofn,jdofn,jnode) 
              
              ! adding elmat because this line computes the contribution of each gauss point
              elmat(ievat,jevat) = elmat(ievat,jevat) &
                   + dvolu * (cn_left * aumat + timas_matrix(idofn,jdofn,jnode,DT_PHYSICAL))


!        esto lo hace mas parecido al explicito antiguo
!              elmat(ievat,jevat) = elmat(ievat,jevat) &
!                   + dvolu * (cn_left * aumat)

              !
              ! In the case of newton or jacobi, timas_matrix(...., DT_PSEUDO) is for the current time step
              ! When tau-newton or tau-jacobi, it is for the real current iter and we have to add the DT_PHYSICAL.
              !
              elmat(ievat,jevat) = elmat(ievat,jevat) &
                   + dvolu * timas_matrix(idofn,jdofn,jnode,DT_PSEUDO)

              
           end if
           
           
           elunk_value = elunk(jdofn,jnode,TIME_N)           
           if (kfl_pseud_nsa == 1) then
              elunk_value = cn_left * elunk(jdofn,jnode,ITER_NEWTON) + cn_pseudo * elunk(jdofn,jnode,TIME_N)
           end if

           if (kfl_linea_nsa==1) then
              ! this is the right hand side contribution K_i * U^n for the delta form Jacobi 
              galte_n(idofn)= galte_n(idofn) &
                   + (advec_matrix(idofn,jdofn,jnode) &
                   +  diffu_matrix(idofn,jdofn,jnode)) * elunk_value
              state_n(idofn)= state_n(idofn) &
                   + (stabi_matrix(idofn,jdofn,jnode)) * elunk_value
              shote_n(idofn)= shote_n(idofn) &
                   + (shote_matrix(idofn,jdofn,jnode)) * elunk_value
              if (kfl_pseud_nsa == 1) then
                 autim_n(idofn) = autim_n(idofn) &
                      + timas_matrix(idofn,jdofn,jnode,DT_PHYSICAL) &
                      * (elunk(jdofn,jnode,ITER_NEWTON) - elunk(jdofn,jnode,TIME_N))
              end if
           else if (kfl_linea_nsa==2) then
              ! this is the right hand side contribution K_i * U^n for the delta form Jacobi 
              autim_n(idofn) = autim_n(idofn) &
                   + timas_matrix(idofn,jdofn,jnode,DT_PSEUDO) &
                   * (elunk(jdofn,jnode,ITER_NEWTON) - elunk(jdofn,jnode,TIME_N))
              galte_n(idofn)= galte_n(idofn) &
                   + (advec_matrix(idofn,jdofn,jnode) &
                   +  diffu_matrix(idofn,jdofn,jnode)) * elunk(jdofn,jnode,ITER_NEWTON)
              state_n(idofn)= state_n(idofn) &
                   + (stabi_matrix(idofn,jdofn,jnode)) * elunk(jdofn,jnode,ITER_NEWTON)
              shote_n(idofn)= shote_n(idofn) &
                   + (shote_matrix(idofn,jdofn,jnode)) * elunk(jdofn,jnode,ITER_NEWTON)
           end if
        end do
     end do

     ! KEEP TIME CONTRIBUTION TO THE RHS TO ZERO WHEN IMPLICIT BECAUSE OF THE DELTA FORM 

     if (kfl_pseud_nsa == 0) autim(idofn) = 0.0_rp

!  esto seria usando xtide
!!!     autim_n(idofn)= xtide(idofn) * xshap(inode)

  end do

  if (kfl_stabi_nsa == 0) then
     state   = 0.0_rp     
     state_n = 0.0_rp     
     shote   = 0.0_rp     
     shote_n = 0.0_rp     
  end if

  if (kfl_timet_nsa == 1) then        
     ! explicit: 
     ! no time contribution to elrhs 
     ! only galte, state and shote explicitly computed terms
     ! elmat is not used
     ! shote comes from outside, it is a positive term
     do idofn = 1,ndaux  
        ievat = (inode-1) * ndofn_nsa + idofn
        ! adding elrhs because this line computes the contribution of each gauss point

!!!!! OJO CON XTIDE ARREGLAAAAAR

        elrhs(ievat) = elrhs(ievat) &
             + dvolu * explicit_pseudo(idofn) * &
             ( - xshap(inode) * xtide(idofn) - galte(idofn) - state(idofn) - shote(idofn))

     end do

     if (ielem==4449 .or. ielem==4500) then
        write(6,*) 'pipit antes',ielem,elrhs             ! el 2
!        write(6,*) advec_matrix(3,1:2,2)
!        write(6,*) 'xc1',xconv(3,1:2,1)
!        write(6,*) 'xc2',xconv(3,1:2,2)
     end if
     if (ielem==4449 .or. ielem==4500) then
        write(6,*) 'pipit despues',ielem,elrhs             ! el 2
!        write(6,*) advec_matrix(3,1:2,2)
!        write(6,*) 'xc1',xconv(3,1:2,1)
!        write(6,*) 'xc2',xconv(3,1:2,2)
     end if

  else if (kfl_timet_nsa == 2) then   
     ! implicit: 
     ! time contribution in elrhs (time n) and elmat (time n+1) 
     ! no galte, state and shote explicitly computed terms
     ! elmat used

!!!     do idofn = 1,ndaux  
!!!        ievat = (inode-1) * ndofn_nsa + idofn
!!!        ! adding elrhs because this line computes the contribution of each gauss point
!!!        elrhs(ievat) = elrhs(ievat) + dvolu * autim(idofn)
!!!     end do

     if (kfl_linea_nsa >= 1) then ! jacobi or inexact newton
        
!        iauxi= 0
!        if (itinn(modul) == 1) then
!           iauxi= 1
!        end if
!        if (iauxi== 1) then
        do idofn = 1,ndaux  
           ievat = (inode-1) * ndofn_nsa + idofn
           resid_n   = 0.0_rp
           if (kfl_linea_nsa == 1) then
              ! jacobi, delta form
              resid_n   = dvolu * (galte_n(idofn) + state_n(idofn) + shote_n(idofn) + autim_n(idofn))
           else if (kfl_linea_nsa == 2) then
              ! newton raphson
              resid_n   = dvolu * (galte_n(idofn) + state_n(idofn) + shote_n(idofn) + autim_n(idofn))
           end if
           elrhs(ievat) = elrhs(ievat) - resid_n
        end do
!        end if
        
     end if
     
  end if
  
!!$
!!$
!!$

!!!!!!!!!!!!!!!!
!!!! OJOOOOO QUE ESTO ES UNA PRUEBA PARA HACER EXPLICITOS!!! NO TIENE QUE IR ASI!!
!!$  auele= 0.0_rp
!!$  do idofn = 1,ndaux  
!!$     ievat = (inode-1) * ndaux + idofn
!!$     do jdofn = 1,ndaux  
!!$        do jnode=1,pnode
!!$           jevat = (jnode-1) * ndofn_nsa + jdofn           
!!$           elunk_jevat  = elunk(jdofn,jnode,ITER_K)
!!$           auele(ievat) = auele(ievat) + elmat(ievat,jevat) * elunk_jevat
!!$        end do
!!$     end do
!!$  end do
!!$  do idofn = 1,ndaux  
!!$     ievat = (inode-1) * ndofn_nsa + idofn
!!$     auele_rhs = elrhs(ievat) + dvolu * ( shote(idofn) - galte(idofn) - state(idofn))
!!$     auele(ievat) = auele(ievat) - auele_rhs
!!$  end do
!!$
!!$  if ((ielem == 2).and.(igaus==4) .and. (inode==4)) then
!!$     write(6,*)
!!$     do ievat=1,nevat_nsa
!!$        write(6,*) auele(ievat)
!!$     end do
!!$     stop
!!$
!!$  end if

!!!!!!!!!!!!!!!!

! para minicucu
!  if (inode == 4 .and. igaus == 4 .and. itinn(modul)== 1) then
!     write(6677,100) itinn(modul),elunk(1,2,ITER_NEWTON), elunk(1,2,TIME_N),autim_n(1),galte_n(1),&
!          cn_left * elunk(1,2,ITER_NEWTON) + cn_pseudo * elunk(1,2,TIME_N)
!     write(6677,100) kfl_pseud_nsa,dtinv_eqs(1,DT_PSEUDO),dtinv_eqs(1,DT_PHYSICAL)
!  end if
!100 format(i,10(2x,e))



end subroutine nsa_elerhs_matrix
