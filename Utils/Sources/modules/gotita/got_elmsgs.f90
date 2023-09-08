subroutine got_elmsgs(&
     pgaus,gpcdr,gpvno,gpdiv,gppor,chale,gpst1,gpst2)
  !------------------------------------------------------------------------
  !****f* Gotita/got_elmsgs
  ! NAME 
  !    got_elmsgs
  ! DESCRIPTION
  !    This routine computes the stabilization parameters 
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_gotita, only       :  kfl_sgsco_got,kfl_sgsti_got,&
       &                        dtinv_got,relsg_got,tosgs_got,&
       &                        itsta_got,resis_got,misgs_got,artif_got,&
       &                        kfl_staty_got,tamin_got,tamax_got,&
       &                        kfl_probl_got,kfl_taust_got,staco_got
  use mod_tauadr, only       :  tauadr
  implicit none
  integer(ip), intent(in)    :: pgaus
  real(rp),    intent(in)    :: chale(2),gpvno(pgaus),gpcdr(pgaus)
  real(rp),    intent(in)    :: gpdiv(pgaus),gppor(pgaus)
  real(rp),    intent(out)   :: gpst1(pgaus),gpst2(pgaus)
  integer(ip)                :: igaus
  real(rp)                   :: freq1,freq2,freq3,freto

  !----------------------------------------------------------------------
  !
  ! BOTH EQUATIONS
  ! 
  !----------------------------------------------------------------------

  if(kfl_probl_got==1) then

     if(kfl_staty_got==1) then         ! SUPG
        do igaus=1,pgaus
           call tauadr(&
                kfl_taust_got,staco_got,gpvno(igaus),0.0_rp,0.0_rp,chale(1),&
                chale(2),gpst1(igaus))
           call tauadr(&
                kfl_taust_got,staco_got,gpvno(igaus),0.0_rp,0.0_rp,chale(1),&
                chale(2),gpst2(igaus))
        end do

     else if(kfl_staty_got>=2) then    ! ASGS and ASGS2
        do igaus=1,pgaus
           call tauadr(&
                kfl_taust_got,staco_got,gpvno(igaus),0.0_rp,gppor(igaus),chale(1),&
                chale(2),gpst1(igaus))
           call tauadr(&
                kfl_taust_got,staco_got,gpvno(igaus),0.0_rp,gpdiv(igaus),chale(1),&
                chale(2),gpst2(igaus))
           !call tauadr(&
           !     kfl_taust_got,staco_got,gpvno(igaus),0.0_rp,0.0_rp,chale(1),&
           !     chale(2),gpst2(igaus))
        end do
     end if

  !----------------------------------------------------------------------
  !
  ! MOMENTUM EQUATION
  ! 
  !----------------------------------------------------------------------

  else if(kfl_probl_got==2) then

     if(kfl_staty_got==1) then         ! SUPG
        do igaus=1,pgaus
           call tauadr(&
                kfl_taust_got,staco_got,gpvno(igaus),0.0_rp,0.0_rp,chale(1),&
                chale(2),gpst1(igaus))
        end do

     else if(kfl_staty_got>=2) then    ! ASGS and ASGS2
        do igaus=1,pgaus
           call tauadr(&
                kfl_taust_got,staco_got,gpvno(igaus),0.0_rp,gppor(igaus),chale(1),&
                chale(2),gpst1(igaus))
        end do

     end if

  !----------------------------------------------------------------------
  !
  ! CONTINUITY EQUATION
  ! 
  !----------------------------------------------------------------------

  else if(kfl_probl_got==3) then

     if(kfl_staty_got==1) then         ! SUPG
        do igaus=1,pgaus
           call tauadr(&
                kfl_taust_got,staco_got,gpvno(igaus),0.0_rp,0.0_rp,chale(1),&
                chale(2),gpst2(igaus))           
        end do

     else if(kfl_staty_got>=2) then    ! ASGS
        do igaus=1,pgaus
           call tauadr(&
                kfl_taust_got,staco_got,gpvno(igaus),0.0_rp,gpdiv(igaus),chale(1),&
                chale(2),gpst2(igaus))           
        end do
     end if

  end if
  !
  ! Minimum and maximum tau
  !
  if(kfl_probl_got==1) then
     do igaus=1,pgaus
        if(gpst1(igaus)>tamax_got) tamax_got=gpst1(igaus)
        if(gpst2(igaus)>tamax_got) tamax_got=gpst2(igaus)
        if(gpst1(igaus)<tamin_got) tamin_got=gpst1(igaus)
        if(gpst2(igaus)<tamin_got) tamin_got=gpst2(igaus)
     end do
  else if(kfl_probl_got==2) then
     do igaus=1,pgaus
        if(gpst1(igaus)>tamax_got) tamax_got=gpst1(igaus)
        if(gpst1(igaus)<tamin_got) tamin_got=gpst1(igaus)
     end do
  else if(kfl_probl_got==3) then
     do igaus=1,pgaus
        if(gpst2(igaus)>tamax_got) tamax_got=gpst2(igaus)
        if(gpst2(igaus)<tamin_got) tamin_got=gpst2(igaus)
     end do
  end if

end subroutine got_elmsgs
