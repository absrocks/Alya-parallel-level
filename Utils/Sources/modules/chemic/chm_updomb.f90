subroutine chm_updomb()
  !------------------------------------------------------------------------
  !****f* partis/chm_updomb
  ! NAME 
  !    chm_updomb
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  implicit none
  integer(ip)       :: ipoin,iodes
  real(rp), pointer :: AA(:,:),CC(:,:),bb(:)
  real(rp)          :: deter
  real(rp)          :: s1,alphamb,betamb,m,alphaw,betaw,Ab,gamma
  real(rp)          :: vf,vw,vl,vf1,vw1,vl1,eta,s2,b,b1
  real(rp), pointer :: denum(:)

  if( nodes_chm /= 0 ) then

     allocate(denum(2*nodes_chm))
     do iodes = 1,nodes_chm
        denum(iodes)           = 0.0_rp ! Numerator
        denum(iodes+nodes_chm) = 0.0_rp ! Denominator
     end do

     if( INOTMASTER ) then

        if( wprob_chm == 'OSTE1' ) then

           !-------------------------------------------------------------------
           !
           ! Problem 1: OSTEOBLASTS 
           !
           ! PDE's: u1 = m, u2 = s1, u3 = s2, u4 = c
           ! ODE's: u5 = b, u6 = vf, u7 = vw, u8 = vl
           !
           !--------------------------------------------------------------------

           alphamb = 0.5_rp    
           betamb  = 10e-6_rp    
           alphaw  = 1.0e-4_rp
           betaw   = 10.0e-6_rp  
           Ab      = 6.67e-3_rp    
           gamma   = 1.0_rp / 365.0_rp

           allocate( AA(3,3) , CC(3,3) , bb(3) )

           do ipoin = 1,npoin

              m   = conce(ipoin,1,1)
              s1  = conce(ipoin,2,1)
              s2  = conce(ipoin,3,1)

              b1  = conce(ipoin,1 + nclas_chm,1)
              vf  = conce(ipoin,2 + nclas_chm,1)
              vw  = conce(ipoin,3 + nclas_chm,1)
              vl  = conce(ipoin,4 + nclas_chm,1)

              b1  = conce(ipoin,1 + nclas_chm,3)
              vf1 = conce(ipoin,2 + nclas_chm,3)
              vw1 = conce(ipoin,3 + nclas_chm,3)
              vl1 = conce(ipoin,4 + nclas_chm,3)

              b       =  ( b1 * dtinv_chm + alphamb * s1 * m / (betamb+s1) ) / ( dtinv_chm + Ab )

              eta     =  b*(alphaw*s2)/(betaw+s2)

              AA(1,1) =  dtinv_chm + eta*(1.0_rp-vw)
              AA(1,2) = -eta * vf
              AA(1,3) =  0.0_rp
              bb(1)   =  dtinv_chm * vf1 - eta*vf*vw

              AA(2,1) = -eta*(1.0_rp-vw) 
              AA(2,2) =  dtinv_chm + eta*vf + gamma*(1.0_rp-vl)
              AA(2,3) = -gamma*vw
              bb(2)   =  dtinv_chm * vw1 + eta*vw*vf - gamma*vw*vl

              AA(3,1) =  0.0_rp
              AA(3,2) = -gamma*(1.0_rp-vl)
              AA(3,3) =  dtinv_chm + gamma*vw
              bb(3)   =  dtinv_chm*vl1 + gamma*vw*vl

              call invmtx(AA,CC,deter,3_ip)

              vf = CC(1,1) * bb(1) + CC(1,2) * bb(2) + CC(1,3) * bb(3)
              vw = CC(2,1) * bb(1) + CC(2,2) * bb(2) + CC(2,3) * bb(3)
              vl = CC(3,1) * bb(1) + CC(3,2) * bb(2) + CC(3,3) * bb(3)
              !
              ! Residual
              !           
              if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then

                 denum(1+nodes_chm) = denum(1+nodes_chm) +  b *  b
                 denum(2+nodes_chm) = denum(2+nodes_chm) + vf * vf
                 denum(3+nodes_chm) = denum(3+nodes_chm) + vw * vw
                 denum(4+nodes_chm) = denum(4+nodes_chm) + vl * vl
                 denum(1)           = denum(1) + (  b - conce(ipoin,1 + nclas_chm,1) )**2
                 denum(2)           = denum(2) + ( vf - conce(ipoin,2 + nclas_chm,1) )**2
                 denum(3)           = denum(3) + ( vw - conce(ipoin,3 + nclas_chm,1) )**2
                 denum(4)           = denum(4) + ( vl - conce(ipoin,4 + nclas_chm,1) )**2

              end if
              !
              ! Solution update
              !           
              conce(ipoin,1 + nclas_chm,1) =  b 
              conce(ipoin,2 + nclas_chm,1) = vf 
              conce(ipoin,3 + nclas_chm,1) = vw 
              conce(ipoin,4 + nclas_chm,1) = vl 

           end do

           deallocate( AA , CC , bb )

        end if

     end if

     call pararr('SUM',0_ip,2*nodes_chm,denum)

     do iodes = 1,nodes_chm
        if( denum(iodes+nodes_chm) /= 0.0_rp ) then
           ripts_chm(iodes + nclas_chm) = sqrt( denum(iodes) / denum(iodes+nodes_chm) )
        end if
        if( INOTSLAVE ) write(lun_spcvg_chm,104,advance='no') ripts_chm(iodes + nclas_chm)
     end do
     if( INOTSLAVE ) then     
        write(lun_spcvg_chm,*) 
        flush(lun_spcvg_chm)
     end if
     deallocate( denum )

  end if

104 format(2x,e12.6)

end subroutine chm_updomb
