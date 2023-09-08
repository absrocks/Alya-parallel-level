subroutine qua_outvar(ivari)
  !------------------------------------------------------------------------
  !****f* Quanty/qua_outvar
  ! NAME 
  !    qua_outvar
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    qua_outvar
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_quanty
  use mod_postpr
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip)             :: ibopo,ipoin,ii
  integer(ip)             :: ieiva
  real(rp)                :: dummr,rutim
  character(5)            :: wopos(3)
  character(2)            :: weiva
  !
  ! Define postprocess variable
  !
  rutim = cutim

  select case (ivari)  

  case(0_ip)
     return

  case(1_ip)
     !
     ! Eigenvalues
     !
     if( INOTSLAVE ) then
        do ieiva=1,eigen_qua(1)%neiva
           !write(lun_witne_qua,*) ieiva,eigva(ieiva)  
        end do
     end if
     !
     ! Eigenvectors
     !
     wopos(2)=wopos_qua(2,1)
     ii = 1
     do ieiva = 1,eigen_qua(1)%neiva
        if( INOTMASTER ) gesca => eigen(ii:) 
        weiva =  intost(ieiva)
        if( ieiva < 10 ) then
           wopos(1) = wopos_qua(1,1)(1:3)//'0'//trim(weiva)
        else
           wopos(1) = wopos_qua(1,1)(1:3)//trim(weiva)
        end if
        call postpr(gesca,wopos,ittim,rutim)
        ii = ii + npoin
     end do
     if( INOTMASTER ) nullify(gesca)
     return

  case(2_ip)

     return
  end select

  ! call outvar(&
  !      ivari,iar3p,&
  !      ittim,rutim,wopos_qua(1,ivari))

end subroutine qua_outvar
