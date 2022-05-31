subroutine skyase(&
     nomat,amatr,nnode,ndime,ndofn,mevat,nequa,&
     lncon,lpdo2,iprob)
!------------------------------------------------------------------------
!    
! This routine performs the assembly for the skyline storage
!
! iprob=1 ... NASTIN
! iprob=2 ... Others
!
! Example for NASTIN. P1/P1 element in 2 dimensions: 2 with 3
!
!       u1 v1   u2 v2   u3 v3   p1 p2 p3
!    +--                               --+
!    |                        .          |
! u1 |  o  o    o  o    o  o  . o  o  o  | 
!    |                        .          |
! v1 |  o  o    o  o    o  o  . o  o  o  |
!    |                        .          |
!    |                 +----+ .      +-+ |
! u2 |  o  o    o  o   |o  o| . o  o |o| |           u3 v3 p3
!    |                 |    | .      | | |          +----+-+
! v2 |  o  o    o  o   |o  o| . o  o |o| |       u2 |o  o|o|
!    |                 +----+ .      +-+ |          |    | |
!    |                        .          |   =>  v2 |o  o|o|
! u3 |  o  o    o  o    o  o  . o  o  o  |          +----+-+
!    |                        .          |       p2 |o  o|o|
! v3 |  o  o    o  o    o  o  . o  o  o  |          +----+-+
!    |                        .          |
!    | ................................. |
! p1 |  o  o    o  o    o  o  . o  o  o  |
!    |                 +----+ .      +-+ |
! p2 |  o  o    o  o   |o  o| . o  o |o| |
!    |                 +----+ .      +-+ |
! p3 |  o  o    o  o    o  o  . o  o  o  |
!    |                        .          |
!    +--                               --+
!
!------------------------------------------------------------------------
  use      def_kintyp
  use      def_solver
  implicit none
  integer(ip), intent(in)    :: nnode,ndime,ndofn,nequa,mevat,iprob
  integer(ip), intent(in)    :: lpdo2(nequa)
  integer(ip), intent(in)    :: lncon(nnode)
  real(rp),    intent(in)    :: nomat(mevat,mevat)
  real(rp),    intent(inout) :: amatr(solve_sol(1)%nzmat)
  integer(ip)                :: idest,jdest,jloca,itot0,jtot0
  integer(ip)                :: itotv,jtotv,idofn,jdofn,i0,i,j,j0
  integer(ip)                :: inode,jnode,idime,jdime,ipoin,jpoin

  select case (iprob)

  case(1)

     do inode=1,nnode
        ipoin=lncon(inode)
        idest=lpntn(ipoin)
        itot0=(idest-1)*ndofn
        i0   =(inode-1)*ndime
        do jnode=1,nnode
           jpoin=lncon(jnode)
           jdest=lpntn(jpoin)
           jtot0=(jdest-1)*ndofn
           j0   =(jnode-1)*ndime
           !
           ! Top left
           !
           do idime = 1,ndime
              itotv = itot0 + idime
              i     = i0    + idime
              do jdime = 1,ndime
                 jtotv = jtot0 + jdime
                 j     = j0    + jdime
                 if(itotv==jtotv) then
                    !$OMP ATOMIC
                    amatr(itotv) = amatr(itotv) + nomat(i,j)
                 else if(itotv<jtotv) then   
                    jloca = itotv + lpdo2(jtotv) - jtotv + 1 + nequa
                    !$OMP ATOMIC
                    amatr(jloca) = amatr(jloca) + nomat(i,j)
                 else 
                    jloca = jtotv + lpdo2(itotv) - itotv + 1 + nequa + lpdo2(nequa)
                    !$OMP ATOMIC
                    amatr(jloca) = amatr(jloca) + nomat(i,j)
                 end if
              end do
           end do

           !
           ! Top right
           !
           itot0 = (idest-1)*ndofn
           jtot0 = (jdest-1)*ndofn
           jtotv = jtot0 + ndime+1
           j     = nnode*ndime+jnode             
           i0    = (inode-1)*ndime

           do idime = 1,ndime
              itotv = itot0 + idime
              i     = i0    + idime
              if(itotv==jtotv) then
                 !$OMP ATOMIC
                 amatr(itotv) = amatr(itotv) + nomat(i,j)
              else if(itotv<jtotv) then   
                 jloca = itotv + lpdo2(jtotv) - jtotv + 1 + nequa
                 !$OMP ATOMIC
                 amatr(jloca) = amatr(jloca) + nomat(i,j)
              else 
                 jloca = jtotv + lpdo2(itotv) - itotv + 1 + nequa + lpdo2(nequa)
                 !$OMP ATOMIC
                 amatr(jloca) = amatr(jloca) + nomat(i,j)
              end if
           end do

           !
           ! Bot. left
           !
           itotv = (idest-1)*ndofn + ndime + 1
           jtot0 = (jdest-1)*ndofn
           i     = nnode*ndime+inode
           j0    = (jnode-1)*ndime       

           do jdime = 1,ndime
              jtotv = jtot0 + jdime
              j     = j0    + jdime        
              if(itotv==jtotv) then
                 !$OMP ATOMIC
                 amatr(itotv) = amatr(itotv) + nomat(i,j)
              else if(itotv<jtotv) then   
                 jloca = itotv + lpdo2(jtotv) - jtotv + 1 + nequa
                 !$OMP ATOMIC
                 amatr(jloca) = amatr(jloca) + nomat(i,j)
              else 
                 jloca = jtotv + lpdo2(itotv) - itotv + 1 + nequa + lpdo2(nequa)
                 !$OMP ATOMIC
                 amatr(jloca) = amatr(jloca) + nomat(i,j)
              end if
           end do

           !
           ! Bot. right
           !
           itotv = (idest-1)*ndofn + ndime + 1
           jtotv = (jdest-1)*ndofn + ndime + 1
           i     = nnode*ndime+inode
           j     = nnode*ndime+jnode            

           if(itotv==jtotv) then
              !$OMP ATOMIC
              amatr(itotv) = amatr(itotv) + nomat(i,j)
           else if(itotv<jtotv) then   
              jloca = itotv + lpdo2(jtotv) - jtotv + 1 + nequa
              !$OMP ATOMIC
              amatr(jloca) = amatr(jloca) + nomat(i,j)
           else 
              jloca = jtotv + lpdo2(itotv) - itotv + 1 + nequa + lpdo2(nequa)
              !$OMP ATOMIC
              amatr(jloca) = amatr(jloca) + nomat(i,j)
           end if

        end do
     end do

  case (2)

     do inode=1,nnode
        itotv=lpntn(lncon(inode))
        do jnode=1,nnode
           jtotv=lpntn(lncon(jnode))

           if(itotv==jtotv) then
              !$OMP ATOMIC
              amatr(itotv) = amatr(itotv) + nomat(inode,jnode)
           else if(itotv<jtotv) then   
              jloca = itotv + lpdo2(jtotv) - jtotv + 1 + nequa
              !$OMP ATOMIC
              amatr(jloca) = amatr(jloca) + nomat(inode,jnode)
           else 
              jloca = jtotv + lpdo2(itotv) - itotv + 1 + nequa + lpdo2(nequa)
              !$OMP ATOMIC
              amatr(jloca) = amatr(jloca) + nomat(inode,jnode)
           end if

        end do
     end do

  case(3)

     do inode=1,nnode
        ipoin = lncon(inode)
        idest = lpntn(ipoin)
        itot0 = (idest-1)*ndofn
        i0    = (inode-1)*ndofn
        do jnode=1,nnode
           jpoin=lncon(jnode)
           jdest=lpntn(jpoin)
           jtot0 = (jdest-1)*ndofn
           j0    = (jnode-1)*ndofn
 
           do idofn=1,ndofn
              itotv=itot0+idofn
              i=i0+idofn
              do jdofn=1,ndofn
                 jtotv = jtot0 + jdofn
                 j=j0+jdofn
                 if(itotv==jtotv) then
                    !$OMP ATOMIC
                    amatr(itotv) = amatr(itotv) + nomat(i,j)
                 else if(itotv<jtotv) then   
                    jloca = itotv + lpdo2(jtotv) - jtotv + 1 + nequa
                    !$OMP ATOMIC
                    amatr(jloca) = amatr(jloca) + nomat(i,j)
                 else 
                    jloca = jtotv + lpdo2(itotv) - itotv + 1 + nequa + lpdo2(nequa)
                    !$OMP ATOMIC
                    amatr(jloca) = amatr(jloca) + nomat(i,j)
                 end if
              end do
           end do

        end do
     end do

  end select

end subroutine skyase

