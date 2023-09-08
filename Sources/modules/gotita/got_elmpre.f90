subroutine got_elmpre(&
     pnode,pgaus,gpsha,gpcar,elcod,elcdr,elvdr,elvel,gpcdr,&
     gpcod,gpvdr,gpvel,gpgvd,gpgcd,gpdiv,gpugu,gpvno)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmpre
  ! NAME 
  !    got_elmpre
  ! DESCRIPTION
  !    Compute some variables at the Gauss points
  ! OUTPUT
  !    GPCOD ... x:           Coordinates
  !    GPVEL ... ua:          Air velocity
  !    GPPRE ... alpha:       Droplet concentration
  !    GPVDR ... u:           Droplet velocity gradients
  !    GPVNO ... |u|:         Droplet Velocity norm
  !    GPGDV ... grad(u):     Droplet velocity gradients
  !    GPDIV ... div(u):      Droplet velocity divergence
  !    GPGCD ... grad(alpha): Droplet concentration gradients
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mgaus,mnode
  use def_gotita, only     :  kfl_timei_got,kfl_probl_got,kfl_linea_got
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: elcod(ndime,pnode),elcdr(pnode,2)
  real(rp),    intent(in)  :: elvdr(ndime,pnode,2),elvel(ndime,pnode)
  real(rp),    intent(out) :: gpcod(ndime,pgaus),gpcdr(pgaus,2)
  real(rp),    intent(out) :: gpvdr(ndime,pgaus,2),gpvel(ndime,pgaus)
  real(rp),    intent(out) :: gpgvd(ndime,ndime,pgaus)
  real(rp),    intent(out) :: gpgcd(ndime,pgaus),gpdiv(pgaus)
  real(rp),    intent(out) :: gpugu(ndime,pgaus),gpvno(pgaus)
  integer(ip)              :: idime,jdime,inode,igaus

  !----------------------------------------------------------------------
  !
  ! BOTH EQUATIONS
  ! 
  !----------------------------------------------------------------------

  if(kfl_probl_got==1) then
     !
     ! Coordinates and air velocity: GPCOD=x, GPVEL=ua
     !
     do igaus=1,pgaus
        do idime=1,ndime
           gpcod(idime,igaus)=0.0_rp
           gpvel(idime,igaus)=0.0_rp
        end do
        do inode=1,pnode
           do idime=1,ndime
              gpcod(idime,igaus)=gpcod(idime,igaus)&
                   +elcod(idime,inode)*gpsha(inode,igaus)
              gpvel(idime,igaus)=gpvel(idime,igaus)&
                   +elvel(idime,inode)*gpsha(inode,igaus)
           end do
        end do
     end do
     !
     ! Water volume fraction and droplet velocity: GPCDR=alpha, GPVDR=u
     !
     do igaus=1,pgaus
        gpcdr(igaus,1)=0.0_rp 
        do idime=1,ndime
           gpvdr(idime,igaus,1)=0.0_rp
        end do
        do inode=1,pnode         
           gpcdr(igaus,1)=gpcdr(igaus,1)&
                +gpsha(inode,igaus)*elcdr(inode,1)
           do idime=1,ndime
              gpvdr(idime,igaus,1)=gpvdr(idime,igaus,1)&
                   +gpsha(inode,igaus)*elvdr(idime,inode,1)
           end do
        end do
     end do
     if(kfl_timei_got==1) then
        do igaus=1,pgaus
           gpcdr(igaus,2)=0.0_rp 
           do idime=1,ndime
              gpvdr(idime,igaus,2)=0.0_rp
           end do
           do inode=1,pnode  
              gpcdr(igaus,2)=gpcdr(igaus,2)&
                   +gpsha(inode,igaus)*elcdr(inode,2)
              do idime=1,ndime
                 gpvdr(idime,igaus,2)=gpvdr(idime,igaus,2)&
                      +gpsha(inode,igaus)*elvdr(idime,inode,2)
              end do
           end do
        end do
     end if
     !
     ! Droplet velocity gradients: GPGVD(j,i)=dui/dxj
     !
     do igaus=1,pgaus
        do idime=1,ndime
           do jdime=1,ndime
              gpgvd(jdime,idime,igaus)=0.0_rp
           end do
        end do
        do inode=1,pnode
           do idime=1,ndime
              do jdime=1,ndime
                 gpgvd(jdime,idime,igaus)=gpgvd(jdime,idime,igaus)&
                      +elvdr(idime,inode,1)*gpcar(jdime,inode,igaus)
              end do
           end do
        end do
     end do
     !
     ! Droplet velocity divergence: GPDIV=div(u)
     !
     do igaus=1,pgaus
        gpdiv(igaus)=0.0_rp
        do idime=1,ndime
           gpdiv(igaus)=gpdiv(igaus)+gpgvd(idime,idime,igaus)
        end do
     end do
     !
     ! Convective term: GPUGU=(u.grad)u
     !
     do igaus=1,pgaus
        do idime=1,ndime
           gpugu(idime,igaus)=0.0_rp
           do jdime=1,ndime
              gpugu(idime,igaus)=gpugu(idime,igaus)&
                   +gpvdr(jdime,igaus,1)*gpgvd(jdime,idime,igaus)
           end do
        end do
     end do
     !
     ! Droplet Velocity norm: GPVNO=|u|
     !
     do igaus=1,pgaus
        call vecnor(gpvdr(1,igaus,1),ndime,gpvno(igaus),2_ip)
     end do
     !
     ! Water volume fraction gradients: GPGCD(i)=d alpha/dxi
     !
     do igaus=1,pgaus
        do idime=1,ndime
           gpgcd(idime,igaus)=0.0_rp  
        end do
     end do
     do inode=1,pnode
        do igaus=1,pgaus
           do idime=1,ndime
              gpgcd(idime,igaus)=gpgcd(idime,igaus)&
                   +gpcar(idime,inode,igaus)*elcdr(inode,1)
           end do
        end do
     end do

  !----------------------------------------------------------------------
  !
  ! DROPLET VELOCITY EQUATION
  ! 
  !----------------------------------------------------------------------

  else if(kfl_probl_got==2) then
     !
     ! Coordinates and air velocity: GPCOD=x, GPVEL=ua
     !
     do igaus=1,pgaus
        do idime=1,ndime
           gpcod(idime,igaus)=0.0_rp
           gpvel(idime,igaus)=0.0_rp
        end do
        do inode=1,pnode
           do idime=1,ndime
              gpcod(idime,igaus)=gpcod(idime,igaus)&
                   +elcod(idime,inode)*gpsha(inode,igaus)
              gpvel(idime,igaus)=gpvel(idime,igaus)&
                   +elvel(idime,inode)*gpsha(inode,igaus)
           end do
        end do
     end do
     !
     ! Water volume fraction and droplet velocity: GPVDR=u
     !
     do igaus=1,pgaus
        do idime=1,ndime
           gpvdr(idime,igaus,1)=0.0_rp
        end do
        do inode=1,pnode
           do idime=1,ndime
              gpvdr(idime,igaus,1)=gpvdr(idime,igaus,1)&
                   +gpsha(inode,igaus)*elvdr(idime,inode,1)
           end do
        end do
     end do
     if(kfl_timei_got==1) then
        do igaus=1,pgaus
           do idime=1,ndime
              gpvdr(idime,igaus,2)=0.0_rp
           end do
           do inode=1,pnode  
              do idime=1,ndime
                 gpvdr(idime,igaus,2)=gpvdr(idime,igaus,2)&
                      +gpsha(inode,igaus)*elvdr(idime,inode,2)
              end do
           end do
        end do
     end if
     !
     ! Droplet velocity gradients: GPGVD(j,i)=dui/dxj
     !
     do igaus=1,pgaus
        do idime=1,ndime
           do jdime=1,ndime
              gpgvd(jdime,idime,igaus)=0.0_rp
           end do
        end do
        do inode=1,pnode
           do idime=1,ndime
              do jdime=1,ndime
                 gpgvd(jdime,idime,igaus)=gpgvd(jdime,idime,igaus)&
                      +elvdr(idime,inode,1)*gpcar(jdime,inode,igaus)
              end do
           end do
        end do
     end do
     !
     ! Droplet velocity divergence: GPDIV=div(u)
     !
     do igaus=1,pgaus
        gpdiv(igaus)=0.0_rp
        do idime=1,ndime
           gpdiv(igaus)=gpdiv(igaus)+gpgvd(idime,idime,igaus)
        end do
     end do
     !
     ! Droplet Velocity norm: GPVNO=|u|
     !
     do igaus=1,pgaus
        call vecnor(gpvdr(1,igaus,1),ndime,gpvno(igaus),2_ip)
     end do
     !
     ! Convective term: GPUGU=(u.grad)u
     !
     if(kfl_linea_got==2) then
        do igaus=1,pgaus
           do idime=1,ndime
              gpugu(idime,igaus)=0.0_rp
              do jdime=1,ndime
                 gpugu(idime,igaus)=gpugu(idime,igaus)&
                      +gpvdr(jdime,igaus,1)*gpgvd(jdime,idime,igaus)
              end do
           end do
        end do
     end if

  !----------------------------------------------------------------------
  !
  ! CONTINUITY EQUATION
  ! 
  !----------------------------------------------------------------------

  else if(kfl_probl_got==3) then
     !
     ! Coordinates and air velocity: GPCOD=x
     !
     do igaus=1,pgaus
        do idime=1,ndime
           gpcod(idime,igaus)=0.0_rp
        end do
        do inode=1,pnode
           do idime=1,ndime
              gpcod(idime,igaus)=gpcod(idime,igaus)&
                   +elcod(idime,inode)*gpsha(inode,igaus)
           end do
        end do
     end do
     !
     ! Water volume fraction and droplet velocity: GPCDR=alpha, GPVDR=u
     !
     do igaus=1,pgaus
        gpcdr(igaus,1)=0.0_rp 
        do idime=1,ndime
           gpvdr(idime,igaus,1)=0.0_rp
        end do
        do inode=1,pnode         
           gpcdr(igaus,1)=gpcdr(igaus,1)&
                +gpsha(inode,igaus)*elcdr(inode,1)
           do idime=1,ndime
              gpvdr(idime,igaus,1)=gpvdr(idime,igaus,1)&
                   +gpsha(inode,igaus)*elvdr(idime,inode,1)
           end do
        end do
     end do
     if(kfl_timei_got==1) then
        do igaus=1,pgaus
           gpcdr(igaus,2)=0.0_rp 
           do inode=1,pnode  
              gpcdr(igaus,2)=gpcdr(igaus,2)&
                   +gpsha(inode,igaus)*elcdr(inode,2)
           end do
        end do
     end if
     !
     ! Droplet velocity gradients: GPGVD(j,i)=dui/dxj
     !
     do igaus=1,pgaus
        do idime=1,ndime
           do jdime=1,ndime
              gpgvd(jdime,idime,igaus)=0.0_rp
           end do
        end do
        do inode=1,pnode
           do idime=1,ndime
              do jdime=1,ndime
                 gpgvd(jdime,idime,igaus)=gpgvd(jdime,idime,igaus)&
                      +elvdr(idime,inode,1)*gpcar(jdime,inode,igaus)
              end do
           end do
        end do
     end do
     !
     ! Droplet velocity divergence: GPDIV=div(u)
     !
     do igaus=1,pgaus
        gpdiv(igaus)=0.0_rp
        do idime=1,ndime
           gpdiv(igaus)=gpdiv(igaus)+gpgvd(idime,idime,igaus)
        end do
     end do
     !
     ! Droplet Velocity norm: GPVNO=|u|
     !
     do igaus=1,pgaus
        call vecnor(gpvdr(1,igaus,1),ndime,gpvno(igaus),2_ip)
     end do

  end if

end subroutine got_elmpre
