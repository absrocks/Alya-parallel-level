subroutine dod_concat(ipoin,isubd,tboun,lext2,nboun_local,ielem_aux)
  use def_parame
  use def_dodeme
  use mod_dod_extens
  use def_master, only : servi,mem_servi
  use def_domain
  use mod_memchk

  implicit none

  integer(ip), intent(in)    :: ipoin,isubd,tboun,nboun_local
  integer(ip), intent(inout) :: ielem_aux,lext2(tboun)
  integer(ip)                :: lpoex_ord(tboun),ii,jj,ipoex,cont,itri,nutri,nucat,ipoiz,imatn
  integer(ip)                :: nfact, index, poex1, poex2,poex3,poex4,flag,tetex,sign
  integer(4)                :: istat
  real(rp)                  :: q_old,minim
  real(rp)   , allocatable  :: q(:)
  integer(ip), allocatable  :: lnods_cat(:,:,:),lext2_tmp(:)
  real(rp)                  :: fact1,fact2,fact3,dumm5(3),kappa,asrad,dummr

  allocate(lext2_tmp(tboun),stat=istat)
  !!???como va??!call memchk(zero,istat,mem_servi(1:2,servi),'lext2_tmp','dod_concat',lext2_tmp)

  !if(tetex/=tboun)then
  !print*,'entro con lista:'
  !print*,lext2(ipoin,1:tboun)
  !print*,'fin-lista'

  !do ii = 1,tboun
  ! lpoex_ord(ii) =lext2(ipoin,ii) 
  ! lext2(ipoin,ii) = 0
  !end do
  !call sortin(tboun,lpoex_ord)
  !cont = 1
  !ipoex = lpoex_ord(1)
  !lext2(ipoin,1)=ipoex
  !do ii=2,tboun
  !  if(lpoex_ord(ii) /= ipoex)thens

  !   cont = cont + 1
  !   lext2(ipoin,cont) = lpoex_ord(ii)
  !   ipoex = lpoex_ord(ii)
  !  end if
  !end do
  ! end if
  !print*,'estoy en el concatttttttt con lista',lext2(ipoin,1:tboun)
  cont = 1
  ipoex=lext2(1)
  lext2_tmp(1) = ipoex
  ii = 0
  do while(ii<=tboun-1)
     ii = ii + 1
     flag = 0
     do jj=1,ii-1
        if(lext2(ii)==lext2(jj))flag=1
     end do
     if(lext2(ii)/=ipoex .and. flag==0)then
        cont = cont + 1
        lext2_tmp(cont)= lext2(ii)
        ipoex = lext2(ii)
     end if
  end do
  flag = 0
  !if(lext2(tboun)/=lext2(tboun-1) .and. lext2(tboun)/=lext2(1))then
  do ii =1,tboun-1
     if(lext2(tboun)==lext2(jj))flag=1  
  end do
  if(flag==0)then
     cont = cont + 1
     lext2_tmp(cont) = lext2(tboun)
     !print*,'he entrado en ultimo if',cont,lext2_tmp(cont)
  end if
  do ii=1,tboun
     lext2(ii) = 0
  end do
  !if(ipoin==147551)print*,'cont',cont
  do ii=1,cont
     lext2(ii)=lext2_tmp(ii)
     ! if(ipoin==147551)print*,lext2(ii)
  end do
  !
  !Numero catalan: numero de posibilidades
  !
  tetex = cont
  if(tetex>2)then
     nutri = tetex - 2
     fact1 = 1
     nfact = 2*nutri
     do ii = 2, nfact
        fact1 = fact1 * ii
     end do

     fact2 = 1
     nfact = nutri
     do ii = 2, nfact
        fact2 = fact2 * ii
     end do

     fact3 = 1
     nfact = nutri+1
     do ii = 2, nfact
        fact3 = fact3 * ii
     end do

     nucat = int(fact1/(fact2*fact3),ip)
  else
     return
  end if

  ! if(DebugMode)print*,'cONCAAAAAAAAAAAAAAAAAAAAAAAAAAT', lext2(1:cont),'ipoinnnn',ipoin,tetex,nucat,nutri

  allocate(q(nucat),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'Q','dod_concat',q)

  allocate(lnods_cat(nucat,nutri,4),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'LNODS_CAT','dod_concat',lnods_cat) 



  if(nucat == 1) then
     ielem_aux=ielem_aux+1
     poex1 = lext2(1)
     poex2 = lext2(2)
     poex3 = lext2(3) 

     !call dod_quaux(ipoin,poex1,poex2,poex3)
     imatn          = lmatn_dod(poex1)
     ipoiz          = lpoiz_dod(poex1)
     lpext(number_fringe_nodes) % ltype(nboun_local+ielem_aux) = 30
     lpext(number_fringe_nodes) % lmate(nboun_local+ielem_aux) = imatn
     lpext(number_fringe_nodes) % lelez(nboun_local+ielem_aux) = ipoiz
     lpext(number_fringe_nodes) % lesub(nboun_local+ielem_aux) = isubd

     call dod_extens_qual3d(3_ip,ipoin,poex1,poex2,poex3,kappa,asrad,dummr,sign,dumm5)
     !if(ipoin== 147551  )
     lpext(number_fringe_nodes) % lnods(1,nboun_local+ielem_aux)= ipoin
     lpext(number_fringe_nodes) % lnods(2,nboun_local+ielem_aux)= poex1
     lpext(number_fringe_nodes) % lnods(3,nboun_local+ielem_aux)= poex2
     lpext(number_fringe_nodes) % lnods(4,nboun_local+ielem_aux)= poex3 
     nelem_dod = nelem_dod + 1
     !print*,'concat1',lpext(number_fringe_nodes) % lnods(:,ielem_aux)

     if(lpext(number_fringe_nodes) % lnods(1,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(2,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(3,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(4,nboun_local+ielem_aux)==0)then
        print*,'STOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOPPPPPP',ielem_aux,ipoin,nucat
        stop
     end if
  else if(nucat==2) then
     ! print*,ipoin,lext2(1),lext2(2),lext2(3),lext2(4)
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(1),sign)   
     q_old = q(1)
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(1),sign)
     if(q(1)<q_old)q(1)=q_old   

     lnods_cat(1,1,1)=ipoin
     lnods_cat(1,1,2)=lext2(1)
     lnods_cat(1,1,3)=lext2(2)
     lnods_cat(1,1,4)=lext2(3)
     lnods_cat(1,2,1)=ipoin
     lnods_cat(1,2,2)=lext2(2)
     lnods_cat(1,2,3)=lext2(3)
     lnods_cat(1,2,4)=lext2(4)


     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(1),q(2),sign)   
     q_old = q(2)
     call dod_qualit(ipoin,lext2(3),lext2(1),lext2(2),q(2),sign)
     if(q(2)<q_old)q(2)=q_old  

     lnods_cat(2,1,1)=ipoin
     lnods_cat(2,1,2)=lext2(3)
     lnods_cat(2,1,3)=lext2(4)
     lnods_cat(2,1,4)=lext2(1)
     lnods_cat(2,2,1)=ipoin
     lnods_cat(2,2,2)=lext2(3)
     lnods_cat(2,2,3)=lext2(1)
     lnods_cat(2,2,4)=lext2(2)

     index = 1
     minim = q(1)

     do ii = 2,nucat
        if(q(ii)<minim ) then
           index = ii
           minim = q(ii)
        end if
     end do

     !if(q(1)<q(2))then
     do itri = 1, nutri
!!!!ielem_aux = -lntyp_dod(ipoin)%l(2*tboun+itri)
        ielem_aux=ielem_aux+1
!!!ltype_aux(ielem_aux) = -104

        poex1 = lnods_cat(index,itri,1)
        poex2 = lnods_cat(index,itri,2)
        poex3 = lnods_cat(index,itri,3)
        poex4 = lnods_cat(index,itri,4)
        imatn          = lmatn_dod(poex1)
        ipoiz          = lpoiz_dod(poex1)
        lpext(number_fringe_nodes) % ltype(nboun_local+ielem_aux) = 30
        lpext(number_fringe_nodes) % lmate(nboun_local+ielem_aux) = imatn
        lpext(number_fringe_nodes) % lelez(nboun_local+ielem_aux) = ipoiz
        lpext(number_fringe_nodes) % lesub(nboun_local+ielem_aux) = isubd
        call dod_extens_qual3d(3_ip,poex1,poex2,poex3,poex4,kappa,asrad,dummr,sign,dumm5)
        !call dod_quaux(poex1,poex2,poex3,poex4)

        lpext(number_fringe_nodes) % lnods(1,nboun_local+ielem_aux)= poex1
        lpext(number_fringe_nodes) % lnods(2,nboun_local+ielem_aux)= poex2
        lpext(number_fringe_nodes) % lnods(3,nboun_local+ielem_aux)= poex3
        lpext(number_fringe_nodes) % lnods(4,nboun_local+ielem_aux)= poex4
        nelem_dod = nelem_dod + 1
        !print*,'concat2',lpext(number_fringe_nodes) % lnods(ielem_aux,:)
        if( lpext(number_fringe_nodes) % lnods(1,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(2,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(3,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(4,nboun_local+ielem_aux)==0 )then
           print*,'STOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOPPPPPP',ielem_aux,ipoin,nucat
           stop
        end if
     end do
  else if(nucat==5) then


     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(1),q(1),sign)   
     q_old = q(1)
     call dod_qualit(ipoin,lext2(4),lext2(1),lext2(3),q(1),sign)
     if(q(1)<q_old)q(1)=q_old   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(1),sign)
     if(q(1)<q_old)q(1)=q_old   

     lnods_cat(1,1,1)=ipoin
     lnods_cat(1,1,2)=lext2(4)
     lnods_cat(1,1,3)=lext2(5)
     lnods_cat(1,1,4)=lext2(1)
     lnods_cat(1,2,1)=ipoin
     lnods_cat(1,2,2)=lext2(4)
     lnods_cat(1,2,3)=lext2(1)
     lnods_cat(1,2,4)=lext2(3)
     lnods_cat(1,3,1)=ipoin
     lnods_cat(1,3,2)=lext2(1)
     lnods_cat(1,3,3)=lext2(2)
     lnods_cat(1,3,4)=lext2(3)


     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(2),sign)   
     q_old = q(2)
     call dod_qualit(ipoin,lext2(3),lext2(5),lext2(2),q(2),sign)
     if(q(2)<q_old)q(2)=q_old  
     call dod_qualit(ipoin,lext2(5),lext2(1),lext2(2),q(2),sign)
     if(q(2)<q_old)q(2)=q_old  

     lnods_cat(2,1,1)=ipoin
     lnods_cat(2,1,2)=lext2(3)
     lnods_cat(2,1,3)=lext2(4)
     lnods_cat(2,1,4)=lext2(5)
     lnods_cat(2,2,1)=ipoin
     lnods_cat(2,2,2)=lext2(3)
     lnods_cat(2,2,3)=lext2(5)
     lnods_cat(2,2,4)=lext2(2)
     lnods_cat(2,3,1)=ipoin
     lnods_cat(2,3,2)=lext2(5)
     lnods_cat(2,3,3)=lext2(1)
     lnods_cat(2,3,4)=lext2(2)
     call dod_qualit(ipoin,lext2(5),lext2(1),lext2(2),q(3),sign)   
     q_old = q(3)
     call dod_qualit(ipoin,lext2(5),lext2(2),lext2(4),q(3),sign)
     if(q(3)<q_old)q(3)=q_old  
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(3),sign)
     if(q(3)<q_old)q(3)=q_old  

     lnods_cat(3,1,1)=ipoin
     lnods_cat(3,1,2)=lext2(5)
     lnods_cat(3,1,3)=lext2(1)
     lnods_cat(3,1,4)=lext2(2)
     lnods_cat(3,2,1)=ipoin
     lnods_cat(3,2,2)=lext2(5)
     lnods_cat(3,2,3)=lext2(2)
     lnods_cat(3,2,4)=lext2(4)
     lnods_cat(3,3,1)=ipoin
     lnods_cat(3,3,2)=lext2(2)
     lnods_cat(3,3,3)=lext2(3)
     lnods_cat(3,3,4)=lext2(4)
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(4),sign)   
     q_old = q(4)
     call dod_qualit(ipoin,lext2(5),lext2(1),lext2(3),q(4),sign)
     if(q(4)<q_old)q(4)=q_old  
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(4),sign)
     if(q(4)<q_old)q(4)=q_old  

     lnods_cat(4,1,1)=ipoin
     lnods_cat(4,1,2)=lext2(3)
     lnods_cat(4,1,3)=lext2(4)
     lnods_cat(4,1,4)=lext2(5)
     lnods_cat(4,2,1)=ipoin
     lnods_cat(4,2,2)=lext2(5)
     lnods_cat(4,2,3)=lext2(1)
     lnods_cat(4,2,4)=lext2(3)
     lnods_cat(4,3,1)=ipoin
     lnods_cat(4,3,2)=lext2(1)
     lnods_cat(4,3,3)=lext2(2)
     lnods_cat(4,3,4)=lext2(3)
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(1),q(5),sign)   
     q_old = q(5)
     call dod_qualit(ipoin,lext2(4),lext2(1),lext2(2),q(5),sign)
     if(q(5)<q_old)q(5)=q_old  
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(5),sign)
     if(q(5)<q_old)q(5)=q_old  

     lnods_cat(5,1,1)=ipoin
     lnods_cat(5,1,2)=lext2(4)
     lnods_cat(5,1,3)=lext2(5)
     lnods_cat(5,1,4)=lext2(1)
     lnods_cat(5,2,1)=ipoin
     lnods_cat(5,2,2)=lext2(4)
     lnods_cat(5,2,3)=lext2(1)
     lnods_cat(5,2,4)=lext2(2)
     lnods_cat(5,3,1)=ipoin
     lnods_cat(5,3,2)=lext2(2)
     lnods_cat(5,3,3)=lext2(3)
     lnods_cat(5,3,4)=lext2(4)

     index = 1
     minim = q(1)

     do ii = 2,nucat
        if(q(ii)<minim ) then
           index = ii
           minim = q(ii)
        end if
     end do
     do itri = 1, nutri
!!! ielem_aux = -lntyp_dod(ipoin)%l(2*tboun+itri)
        !!ltype_ext(ielem_aux) = -104

        ielem_aux             = ielem_aux + 1
        poex1=lnods_cat(index,itri,1)
        poex2=lnods_cat(index,itri,2)
        poex3=lnods_cat(index,itri,3)
        poex4=lnods_cat(index,itri,4)
        imatn          = lmatn_dod(poex1)
        ipoiz          = lpoiz_dod(poex1)
        lpext(number_fringe_nodes) % ltype(nboun_local+ielem_aux) = 30
        lpext(number_fringe_nodes) % lmate(nboun_local+ielem_aux) = imatn
        lpext(number_fringe_nodes) % lelez(nboun_local+ielem_aux) = ipoiz
        lpext(number_fringe_nodes) % lesub(nboun_local+ielem_aux) = isubd
        call dod_extens_qual3d(3_ip,poex1,poex2,poex3,poex4,kappa,asrad,dummr,sign,dumm5)
        !call dod_quaux(poex1,poex2,poex3,poex4)

        lpext(number_fringe_nodes) % lnods(1,nboun_local+ielem_aux)= poex1
        lpext(number_fringe_nodes) % lnods(2,nboun_local+ielem_aux)= poex2
        lpext(number_fringe_nodes) % lnods(3,nboun_local+ielem_aux)= poex3
        lpext(number_fringe_nodes) % lnods(4,nboun_local+ielem_aux)= poex4
        nelem_dod = nelem_dod + 1
        !print*,'concat3',lpext(number_fringe_nodes) % lnods(ielem_aux,:)
        if(lpext(number_fringe_nodes) % lnods(1,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(2,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(3,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(4,nboun_local+ielem_aux)==0)then
           print*,'STOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOPPPPPP',ielem_aux,ipoin,nucat
           stop
        end if
     end do

  else if(nucat==14) then
     print*,'estoy aquiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii',lext2(1:6)
     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(1),q(1),sign)   
     q_old = q(1)
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(1),q(1),sign)
     if(q(1)<q_old)q(1)=q_old   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(4),q(1),sign)
     if(q(1)<q_old)q(1)=q_old   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(1),sign)
     if(q(1)<q_old)q(1)=q_old   


     lnods_cat(1,1,1)=ipoin
     lnods_cat(1,1,2)=lext2(5)
     lnods_cat(1,1,3)=lext2(6)
     lnods_cat(1,1,4)=lext2(1)
     lnods_cat(1,2,1)=ipoin
     lnods_cat(1,2,2)=lext2(4)
     lnods_cat(1,2,3)=lext2(5)
     lnods_cat(1,2,4)=lext2(1)
     lnods_cat(1,3,1)=ipoin
     lnods_cat(1,3,2)=lext2(1)
     lnods_cat(1,3,3)=lext2(2)
     lnods_cat(1,3,4)=lext2(4)
     lnods_cat(1,4,1)=ipoin
     lnods_cat(1,4,2)=lext2(2)
     lnods_cat(1,4,3)=lext2(3)
     lnods_cat(1,4,4)=lext2(4)


     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(2),sign)   
     q_old = q(2)
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(6),q(2),sign)
     if(q(2)<q_old)q(2)=q_old  
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(4),q(2),sign)
     if(q(2)<q_old)q(2)=q_old  
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(2),sign)
     if(q(2)<q_old)q(2)=q_old  

     lnods_cat(2,1,1)=ipoin
     lnods_cat(2,1,2)=lext2(4)
     lnods_cat(2,1,3)=lext2(5)
     lnods_cat(2,1,4)=lext2(6)
     lnods_cat(2,2,1)=ipoin
     lnods_cat(2,2,2)=lext2(1)
     lnods_cat(2,2,3)=lext2(4)
     lnods_cat(2,2,4)=lext2(6)
     lnods_cat(2,3,1)=ipoin
     lnods_cat(2,3,2)=lext2(1)
     lnods_cat(2,3,3)=lext2(2)
     lnods_cat(2,3,4)=lext2(4)
     lnods_cat(2,4,1)=ipoin
     lnods_cat(2,4,2)=lext2(1)
     lnods_cat(2,4,3)=lext2(2)
     lnods_cat(2,4,4)=lext2(4)

     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(6),q(3),sign)   
     q_old = q(3)
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(5),q(3),sign)
     if(q(3)<q_old)q(3)=q_old  
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(5),q(3),sign)
     if(q(3)<q_old)q(3)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(3),sign)
     if(q(3)<q_old)q(3)=q_old  

     lnods_cat(3,1,1)=ipoin
     lnods_cat(3,1,2)=lext2(1)
     lnods_cat(3,1,3)=lext2(5)
     lnods_cat(3,1,4)=lext2(6)
     lnods_cat(3,2,1)=ipoin
     lnods_cat(3,2,2)=lext2(1)
     lnods_cat(3,2,3)=lext2(2)
     lnods_cat(3,2,4)=lext2(5)
     lnods_cat(3,3,1)=ipoin
     lnods_cat(3,3,2)=lext2(2)
     lnods_cat(3,3,3)=lext2(3)
     lnods_cat(3,3,4)=lext2(5)
     lnods_cat(3,4,1)=ipoin
     lnods_cat(3,4,2)=lext2(3)
     lnods_cat(3,4,3)=lext2(4)
     lnods_cat(3,4,4)=lext2(5)

     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(6),q(4),sign)   
     q_old = q(4)
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(6),q(4),sign)
     if(q(4)<q_old)q(4)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(6),q(4),sign)
     if(q(4)<q_old)q(4)=q_old  
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(4),sign)
     if(q(4)<q_old)q(4)=q_old  

     lnods_cat(4,1,1)=ipoin
     lnods_cat(4,1,2)=lext2(1)
     lnods_cat(4,1,3)=lext2(2)
     lnods_cat(4,1,4)=lext2(6)
     lnods_cat(4,2,1)=ipoin
     lnods_cat(4,2,2)=lext2(2)
     lnods_cat(4,2,3)=lext2(3)
     lnods_cat(4,2,4)=lext2(6)
     lnods_cat(4,3,1)=ipoin
     lnods_cat(4,3,2)=lext2(3)
     lnods_cat(4,3,3)=lext2(4)
     lnods_cat(4,3,4)=lext2(6)
     lnods_cat(4,4,1)=ipoin
     lnods_cat(4,4,2)=lext2(4)
     lnods_cat(4,4,3)=lext2(5)
     lnods_cat(4,4,4)=lext2(6)

     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(6),q(5),sign)   
     q_old = q(5)
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(5),q(5),sign)
     if(q(5)<q_old)q(5)=q_old  
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(4),q(5),sign)
     if(q(5)<q_old)q(5)=q_old  
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(5),sign)
     if(q(5)<q_old)q(5)=q_old  

     lnods_cat(5,1,1)=ipoin
     lnods_cat(5,1,2)=lext2(1)
     lnods_cat(5,1,3)=lext2(5)
     lnods_cat(5,1,4)=lext2(6)
     lnods_cat(5,2,1)=ipoin
     lnods_cat(5,2,2)=lext2(1)
     lnods_cat(5,2,3)=lext2(4)
     lnods_cat(5,2,4)=lext2(5)
     lnods_cat(5,3,1)=ipoin
     lnods_cat(5,3,2)=lext2(1)
     lnods_cat(5,3,3)=lext2(3)
     lnods_cat(5,3,4)=lext2(4)
     lnods_cat(5,4,1)=ipoin
     lnods_cat(5,4,2)=lext2(1)
     lnods_cat(5,4,3)=lext2(2)
     lnods_cat(5,4,4)=lext2(3)


     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(6),q(6),sign)   
     q_old = q(6)
     call dod_qualit(ipoin,lext2(2),lext2(5),lext2(6),q(6),sign)
     if(q(6)<q_old)q(6)=q_old   
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(5),q(6),sign)
     if(q(6)<q_old)q(6)=q_old   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(6),sign)
     if(q(6)<q_old)q(6)=q_old   

     lnods_cat(6,1,1)=ipoin
     lnods_cat(6,1,2)=lext2(1)
     lnods_cat(6,1,3)=lext2(2)
     lnods_cat(6,1,4)=lext2(6)
     lnods_cat(6,2,1)=ipoin
     lnods_cat(6,2,2)=lext2(2)
     lnods_cat(6,2,3)=lext2(5)
     lnods_cat(6,2,4)=lext2(6)
     lnods_cat(6,3,1)=ipoin
     lnods_cat(6,3,2)=lext2(2)
     lnods_cat(6,3,3)=lext2(4)
     lnods_cat(6,3,4)=lext2(5)
     lnods_cat(6,4,1)=ipoin
     lnods_cat(6,4,2)=lext2(2)
     lnods_cat(6,4,3)=lext2(3)
     lnods_cat(6,4,4)=lext2(4)


     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(7),sign)   
     q_old = q(7)
     call dod_qualit(ipoin,lext2(3),lext2(5),lext2(6),q(7),sign)
     if( q(7) < q_old ) q(7) = q_old  
     call dod_qualit(ipoin,lext2(3),lext2(6),lext2(1),q(7),sign)
     if( q(7) < q_old ) q(7) = q_old  
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(7),sign)
     if( q(7) < q_old ) q(7) = q_old  

     lnods_cat(7,1,1)=ipoin
     lnods_cat(7,1,2)=lext2(3)
     lnods_cat(7,1,3)=lext2(4)
     lnods_cat(7,1,4)=lext2(5)
     lnods_cat(7,2,1)=ipoin
     lnods_cat(7,2,2)=lext2(3)
     lnods_cat(7,2,3)=lext2(5)
     lnods_cat(7,2,4)=lext2(6)
     lnods_cat(7,3,1)=ipoin
     lnods_cat(7,3,2)=lext2(3)
     lnods_cat(7,3,3)=lext2(6)
     lnods_cat(7,3,4)=lext2(1)
     lnods_cat(7,4,1)=ipoin
     lnods_cat(7,4,2)=lext2(1)
     lnods_cat(7,4,3)=lext2(2)
     lnods_cat(7,4,4)=lext2(3)

     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(6),q(8),sign)   
     q_old = q(8)
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(5),q(8),sign)
     if(q(8)<q_old)q(8)=q_old  
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(8),sign)
     if(q(8)<q_old)q(8)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(8),sign)
     if(q(8)<q_old)q(8)=q_old  

     lnods_cat(8,1,1)=ipoin
     lnods_cat(8,1,2)=lext2(1)
     lnods_cat(8,1,3)=lext2(5)
     lnods_cat(8,1,4)=lext2(6)
     lnods_cat(8,2,1)=ipoin
     lnods_cat(8,2,2)=lext2(1)
     lnods_cat(8,2,3)=lext2(3)
     lnods_cat(8,2,4)=lext2(5)
     lnods_cat(8,3,1)=ipoin
     lnods_cat(8,3,2)=lext2(1)
     lnods_cat(8,3,3)=lext2(2)
     lnods_cat(8,3,4)=lext2(3)
     lnods_cat(8,4,1)=ipoin
     lnods_cat(8,4,2)=lext2(3)
     lnods_cat(8,4,3)=lext2(4)
     lnods_cat(8,4,4)=lext2(5)

     call dod_qualit(ipoin,lext2(6),lext2(1),lext2(2),q(9),sign)   
     q_old = q(9)
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(9),sign)
     if(q(9)<q_old)q(9)=q_old  
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(6),q(9),sign)
     if(q(9)<q_old)q(9)=q_old  
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(9),sign)
     if(q(9)<q_old)q(9)=q_old  

     lnods_cat(9,1,1)=ipoin
     lnods_cat(9,1,2)=lext2(6)
     lnods_cat(9,1,3)=lext2(1)
     lnods_cat(9,1,4)=lext2(2)
     lnods_cat(9,2,1)=ipoin
     lnods_cat(9,2,2)=lext2(4)
     lnods_cat(9,2,3)=lext2(5)
     lnods_cat(9,2,4)=lext2(6)
     lnods_cat(9,3,1)=ipoin
     lnods_cat(9,3,2)=lext2(2)
     lnods_cat(9,3,3)=lext2(4)
     lnods_cat(9,3,4)=lext2(6)
     lnods_cat(9,4,1)=ipoin
     lnods_cat(9,4,2)=lext2(2)
     lnods_cat(9,4,3)=lext2(3)
     lnods_cat(9,4,4)=lext2(4)

     call dod_qualit(ipoin,lext2(6),lext2(1),lext2(2),q(10),sign)   
     q_old = q(10)
     call dod_qualit(ipoin,lext2(6),lext2(2),lext2(3),q(10),sign)
     if(q(10)<q_old)q(10)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(5),lext2(6),q(10),sign)
     if(q(10)<q_old)q(10)=q_old  
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(3),q(10),sign)
     if(q(10)<q_old)q(10)=q_old  

     lnods_cat(10,1,1)=ipoin
     lnods_cat(10,1,2)=lext2(6)
     lnods_cat(10,1,3)=lext2(1)
     lnods_cat(10,1,4)=lext2(2)
     lnods_cat(10,2,1)=ipoin
     lnods_cat(10,2,2)=lext2(6)
     lnods_cat(10,2,3)=lext2(2)
     lnods_cat(10,2,4)=lext2(3)
     lnods_cat(10,3,1)=ipoin
     lnods_cat(10,3,2)=lext2(3)
     lnods_cat(10,3,3)=lext2(4)
     lnods_cat(10,3,4)=lext2(5)
     lnods_cat(10,4,1)=ipoin
     lnods_cat(10,4,2)=lext2(4)
     lnods_cat(10,4,3)=lext2(5)
     lnods_cat(10,4,4)=lext2(3)

     call dod_qualit(ipoin,lext2(6),lext2(1),lext2(2),q(11),sign)   
     q_old = q(11)
     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(2),q(11),sign)
     if(q(11)<q_old)q(11)=q_old  
     call dod_qualit(ipoin,lext2(5),lext2(2),lext2(3),q(11),sign)
     if(q(11)<q_old)q(11)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(11),sign)
     if(q(11)<q_old)q(11)=q_old  

     lnods_cat(11,1,1)=ipoin
     lnods_cat(11,1,2)=lext2(6)
     lnods_cat(11,1,3)=lext2(1)
     lnods_cat(11,1,4)=lext2(2)
     lnods_cat(11,2,1)=ipoin
     lnods_cat(11,2,2)=lext2(5)
     lnods_cat(11,2,3)=lext2(6)
     lnods_cat(11,2,4)=lext2(2)
     lnods_cat(11,3,1)=ipoin
     lnods_cat(11,3,2)=lext2(5)
     lnods_cat(11,3,3)=lext2(2)
     lnods_cat(11,3,4)=lext2(3)
     lnods_cat(11,4,1)=ipoin
     lnods_cat(11,4,2)=lext2(3)
     lnods_cat(11,4,3)=lext2(4)
     lnods_cat(11,4,4)=lext2(5)

     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(12),sign)   
     q_old = q(12)
     call dod_qualit(ipoin,lext2(4),lext2(6),lext2(3),q(12),sign)
     if(q(12)<q_old)q(12)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(6),lext2(1),q(12),sign)
     if(q(12)<q_old)q(12)=q_old  
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(12),sign)
     if(q(12)<q_old)q(12)=q_old  

     lnods_cat(12,1,1)=ipoin
     lnods_cat(12,1,2)=lext2(4)
     lnods_cat(12,1,3)=lext2(5)
     lnods_cat(12,1,4)=lext2(6)
     lnods_cat(12,2,1)=ipoin
     lnods_cat(12,2,2)=lext2(4)
     lnods_cat(12,2,3)=lext2(6)
     lnods_cat(12,2,4)=lext2(3)
     lnods_cat(12,3,1)=ipoin
     lnods_cat(12,3,2)=lext2(3)
     lnods_cat(12,3,3)=lext2(6)
     lnods_cat(12,3,4)=lext2(1)
     lnods_cat(12,4,1)=ipoin
     lnods_cat(12,4,2)=lext2(1)
     lnods_cat(12,4,3)=lext2(2)
     lnods_cat(12,4,4)=lext2(3)


     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(13),sign)   
     q_old = q(13)
     call dod_qualit(ipoin,lext2(4),lext2(6),lext2(1),q(13),sign)
     if(q(13)<q_old)q(13)=q_old  
     call dod_qualit(ipoin,lext2(4),lext2(1),lext2(3),q(13),sign)
     if(q(13)<q_old)q(13)=q_old  
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(13),sign)
     if(q(13)<q_old)q(13)=q_old  

     lnods_cat(13,1,1)=ipoin
     lnods_cat(13,1,2)=lext2(4)
     lnods_cat(13,1,3)=lext2(5)
     lnods_cat(13,1,4)=lext2(6)
     lnods_cat(13,2,1)=ipoin
     lnods_cat(13,2,2)=lext2(4)
     lnods_cat(13,2,3)=lext2(6)
     lnods_cat(13,2,4)=lext2(1)
     lnods_cat(13,3,1)=ipoin
     lnods_cat(13,3,2)=lext2(4)
     lnods_cat(13,3,3)=lext2(1)
     lnods_cat(13,3,4)=lext2(3)
     lnods_cat(13,4,1)=ipoin
     lnods_cat(13,4,2)=lext2(1)
     lnods_cat(13,4,3)=lext2(2)
     lnods_cat(13,4,4)=lext2(3)
     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(6),q(14),sign)   
     q_old = q(14)
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(5),q(14),sign)
     if(q(14)<q_old)q(14)=q_old  
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(5),q(14),sign)
     if(q(14)<q_old)q(14)=q_old  
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(14),sign)
     if(q(14)<q_old)q(14)=q_old  

     lnods_cat(14,1,1)=ipoin
     lnods_cat(14,1,2)=lext2(1)
     lnods_cat(14,1,3)=lext2(5)
     lnods_cat(14,1,4)=lext2(6)
     lnods_cat(14,2,1)=ipoin
     lnods_cat(14,2,2)=lext2(1)
     lnods_cat(14,2,3)=lext2(2)
     lnods_cat(14,2,4)=lext2(5)
     lnods_cat(14,3,1)=ipoin
     lnods_cat(14,3,2)=lext2(2)
     lnods_cat(14,3,3)=lext2(4)
     lnods_cat(14,3,4)=lext2(5)
     lnods_cat(14,4,1)=ipoin
     lnods_cat(14,4,2)=lext2(2)
     lnods_cat(14,4,3)=lext2(3)
     lnods_cat(14,4,4)=lext2(4)

     index = 1
     minim = q(1)

     do ii = 2,nucat
        if(q(ii)<minim ) then
           index = ii
           minim = q(ii)
        end if
     end do
     do itri = 1, nutri
!!!OJOOOO ielem_aux = -lntyp_dod(ipoin)%l(2*tboun+itri)
        ielem_aux             = ielem_aux + 1
!!!ltype_ext(ielem_aux) = -104
        poex1=lnods_cat(index,itri,1)
        poex2=lnods_cat(index,itri,2)
        poex3=lnods_cat(index,itri,3)
        poex4=lnods_cat(index,itri,4)
        imatn          = lmatn_dod(poex1)
        ipoiz          = lpoiz_dod(poex1)
        lpext(number_fringe_nodes) % ltype(nboun_local+ielem_aux) = 30
        lpext(number_fringe_nodes) % lmate(nboun_local+ielem_aux) = imatn
        lpext(number_fringe_nodes) % lelez(nboun_local+ielem_aux) = ipoiz
        lpext(number_fringe_nodes) % lesub(nboun_local+ielem_aux) = isubd
        call dod_extens_qual3d(3_ip,poex1,poex2,poex3,poex4,kappa,asrad,dummr,sign,dumm5)
        !call dod_quaux(poex1,poex2,poex3,poex4)

        lpext(number_fringe_nodes) % lnods(1,nboun_local+ielem_aux)= poex1
        lpext(number_fringe_nodes) % lnods(2,nboun_local+ielem_aux)= poex2
        lpext(number_fringe_nodes) % lnods(3,nboun_local+ielem_aux)= poex3
        lpext(number_fringe_nodes) % lnods(4,nboun_local+ielem_aux)= poex4
        nelem_dod = nelem_dod + 1
        !print*,'concat4',lpext(number_fringe_nodes) % lnods(ielem_aux,:)
        ! if(ipoin==15965)print*,lpext(number_fringe_nodes) % lnods(ielem_aux,:)
        if(lpext(number_fringe_nodes) % lnods(1,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(2,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(3,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(4,nboun_local+ielem_aux)==0)then
           print*,'STOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOPPPPPP',ielem_aux,ipoin,nucat
           stop
        end if
     end do


  else if(nucat==42) then 

     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(1),sign)   
     q_old = q(1)
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(7),q(1),sign)
     if(q(1)<q_old)q(1)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(5),lext2(6),q(1),sign)
     if(q(1)<q_old)q(1)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(5),lext2(6),q(1),sign)
     if(q(1)<q_old)q(1)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(1),sign)
     if(q(1)<q_old)q(1)=q_old  

     lnods_cat(1,1,1)=ipoin
     lnods_cat(1,1,2)=lext2(1)
     lnods_cat(1,1,3)=lext2(2)
     lnods_cat(1,1,4)=lext2(3)
     lnods_cat(1,2,1)=ipoin
     lnods_cat(1,2,2)=lext2(1)
     lnods_cat(1,2,3)=lext2(3)
     lnods_cat(1,2,4)=lext2(7)
     lnods_cat(1,3,1)=ipoin
     lnods_cat(1,3,2)=lext2(3)
     lnods_cat(1,3,3)=lext2(6)
     lnods_cat(1,3,4)=lext2(7)
     lnods_cat(1,4,1)=ipoin
     lnods_cat(1,4,2)=lext2(3)
     lnods_cat(1,4,3)=lext2(5)
     lnods_cat(1,4,4)=lext2(6)
     lnods_cat(1,5,1)=ipoin
     lnods_cat(1,5,2)=lext2(3)
     lnods_cat(1,5,3)=lext2(4)
     lnods_cat(1,5,4)=lext2(5)

     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(2),sign)   
     q_old = q(2)
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(7),q(2),sign)
     if(q(2)<q_old)q(2)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(6),lext2(7),q(2),sign)
     if(q(2)<q_old)q(2)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(5),lext2(6),q(2),sign)
     if(q(2)<q_old)q(2)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(2),sign)
     if(q(2)<q_old)q(2)=q_old  

     lnods_cat(2,1,1)=ipoin
     lnods_cat(2,1,2)=lext2(1)
     lnods_cat(2,1,3)=lext2(2)
     lnods_cat(2,1,4)=lext2(7)
     lnods_cat(2,2,1)=ipoin
     lnods_cat(2,2,2)=lext2(2)
     lnods_cat(2,2,3)=lext2(3)
     lnods_cat(2,2,4)=lext2(7)
     lnods_cat(2,3,1)=ipoin
     lnods_cat(2,3,2)=lext2(3)
     lnods_cat(2,3,3)=lext2(6)
     lnods_cat(2,3,4)=lext2(7)
     lnods_cat(2,4,1)=ipoin
     lnods_cat(2,4,2)=lext2(3)
     lnods_cat(2,4,3)=lext2(5)
     lnods_cat(2,4,4)=lext2(6)
     lnods_cat(2,5,1)=ipoin
     lnods_cat(2,5,2)=lext2(3)
     lnods_cat(2,5,3)=lext2(4)
     lnods_cat(2,5,4)=lext2(5)

     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(3),sign)   
     q_old = q(3)
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(6),q(3),sign)
     if(q(3)<q_old)q(3)=q_old  
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(3),sign)
     if(q(3)<q_old)q(3)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(5),lext2(6),q(3),sign)
     if(q(3)<q_old)q(3)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(3),sign)
     if(q(3)<q_old)q(3)=q_old  
     lnods_cat(3,1,1)=ipoin
     lnods_cat(3,1,2)=lext2(1)
     lnods_cat(3,1,3)=lext2(6)
     lnods_cat(3,1,4)=lext2(7)
     lnods_cat(3,2,1)=ipoin
     lnods_cat(3,2,2)=lext2(1)
     lnods_cat(3,2,3)=lext2(3)
     lnods_cat(3,2,4)=lext2(6)
     lnods_cat(3,3,1)=ipoin
     lnods_cat(3,3,2)=lext2(1)
     lnods_cat(3,3,3)=lext2(2)
     lnods_cat(3,3,4)=lext2(3)
     lnods_cat(3,4,1)=ipoin
     lnods_cat(3,4,2)=lext2(3)
     lnods_cat(3,4,3)=lext2(5)
     lnods_cat(3,4,4)=lext2(6)
     lnods_cat(3,5,1)=ipoin
     lnods_cat(3,5,2)=lext2(3)
     lnods_cat(3,5,3)=lext2(4)
     lnods_cat(3,5,4)=lext2(5)


     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(4),sign)   
     q_old = q(4)
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(6),q(4),sign)
     if(q(4)<q_old)q(4)=q_old  
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(6),q(4),sign)
     if(q(4)<q_old)q(4)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(5),lext2(6),q(4),sign)
     if(q(4)<q_old)q(4)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(4),sign)
     if(q(4)<q_old)q(4)=q_old  
     lnods_cat(4,1,1)=ipoin
     lnods_cat(4,1,2)=lext2(1)
     lnods_cat(4,1,3)=lext2(6)
     lnods_cat(4,1,4)=lext2(7)
     lnods_cat(4,2,1)=ipoin
     lnods_cat(4,2,2)=lext2(1)
     lnods_cat(4,2,3)=lext2(2)
     lnods_cat(4,2,4)=lext2(6)
     lnods_cat(4,3,1)=ipoin
     lnods_cat(4,3,2)=lext2(2)
     lnods_cat(4,3,3)=lext2(3)
     lnods_cat(4,3,4)=lext2(6)
     lnods_cat(4,4,1)=ipoin
     lnods_cat(4,4,2)=lext2(3)
     lnods_cat(4,4,3)=lext2(5)
     lnods_cat(4,4,4)=lext2(6)
     lnods_cat(4,5,1)=ipoin
     lnods_cat(4,5,2)=lext2(3)
     lnods_cat(4,5,3)=lext2(4)
     lnods_cat(4,5,4)=lext2(5)


     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(5),sign)   
     q_old = q(5)
     call dod_qualit(ipoin,lext2(2),lext2(6),lext2(7),q(5),sign)
     if(q(5)<q_old)q(5)=q_old  
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(6),q(5),sign)
     if(q(5)<q_old)q(5)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(5),lext2(6),q(5),sign)
     if(q(5)<q_old)q(5)=q_old  
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(5),sign)
     if(q(5)<q_old)q(5)=q_old  
     lnods_cat(5,1,1)=ipoin
     lnods_cat(5,1,2)=lext2(1)
     lnods_cat(5,1,3)=lext2(2)
     lnods_cat(5,1,4)=lext2(7)
     lnods_cat(5,2,1)=ipoin
     lnods_cat(5,2,2)=lext2(2)
     lnods_cat(5,2,3)=lext2(6)
     lnods_cat(5,2,4)=lext2(7)
     lnods_cat(5,3,1)=ipoin
     lnods_cat(5,3,2)=lext2(2)
     lnods_cat(5,3,3)=lext2(3)
     lnods_cat(5,3,4)=lext2(6)
     lnods_cat(5,4,1)=ipoin
     lnods_cat(5,4,2)=lext2(3)
     lnods_cat(5,4,3)=lext2(5)
     lnods_cat(5,4,4)=lext2(6)
     lnods_cat(5,5,1)=ipoin
     lnods_cat(5,5,2)=lext2(3)
     lnods_cat(5,5,3)=lext2(4)
     lnods_cat(5,5,4)=lext2(5)


     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(6),sign)   
     q_old = q(6)                                                               
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(6),sign)
     if(q(6)<q_old)q(6)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(5),lext2(7),q(6),sign)
     if(q(6)<q_old)q(6)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(7),q(6),sign)
     if(q(6)<q_old)q(6)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(6),sign)
     if(q(6)<q_old)q(6)=q_old  
     lnods_cat(6,1,1)=ipoin
     lnods_cat(6,1,2)=lext2(5)
     lnods_cat(6,1,3)=lext2(6)
     lnods_cat(6,1,4)=lext2(7)
     lnods_cat(6,2,1)=ipoin
     lnods_cat(6,2,2)=lext2(3)
     lnods_cat(6,2,3)=lext2(4)
     lnods_cat(6,2,4)=lext2(5)
     lnods_cat(6,3,1)=ipoin
     lnods_cat(6,3,2)=lext2(3)
     lnods_cat(6,3,3)=lext2(5)
     lnods_cat(6,3,4)=lext2(7)
     lnods_cat(6,4,1)=ipoin
     lnods_cat(6,4,2)=lext2(1)
     lnods_cat(6,4,3)=lext2(3)
     lnods_cat(6,4,4)=lext2(7)
     lnods_cat(6,5,1)=ipoin
     lnods_cat(6,5,2)=lext2(1)
     lnods_cat(6,5,3)=lext2(2)
     lnods_cat(6,5,4)=lext2(3)

     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(7),sign)   
     q_old = q(7)                                                                
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(7),q(7),sign)
     if(q(7)<q_old)q(7)=q_old                                                    
     call dod_qualit(ipoin,lext2(3),lext2(5),lext2(7),q(7),sign)
     if(q(7)<q_old)q(7)=q_old                                                    
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(7),sign)
     if(q(7)<q_old)q(7)=q_old                                                    
     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(7),sign)
     if(q(7)<q_old)q(7)=q_old  
     lnods_cat(7,1,1)=ipoin
     lnods_cat(7,1,2)=lext2(1)
     lnods_cat(7,1,3)=lext2(2)
     lnods_cat(7,1,4)=lext2(7)
     lnods_cat(7,2,1)=ipoin
     lnods_cat(7,2,2)=lext2(2)
     lnods_cat(7,2,3)=lext2(3)
     lnods_cat(7,2,4)=lext2(7)
     lnods_cat(7,3,1)=ipoin
     lnods_cat(7,3,2)=lext2(3)
     lnods_cat(7,3,3)=lext2(5)
     lnods_cat(7,3,4)=lext2(7)
     lnods_cat(7,4,1)=ipoin
     lnods_cat(7,4,2)=lext2(3)
     lnods_cat(7,4,3)=lext2(4)
     lnods_cat(7,4,4)=lext2(5)
     lnods_cat(7,5,1)=ipoin
     lnods_cat(7,5,2)=lext2(5)
     lnods_cat(7,5,3)=lext2(6)
     lnods_cat(7,5,4)=lext2(7)


     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(8),sign)   
     q_old = q(8)                                                                
     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(7),q(8),sign)
     if(q(8)<q_old)q(8)=q_old                                                    
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(5),q(8),sign)
     if(q(8)<q_old)q(8)=q_old                                                    
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(8),sign)
     if(q(8)<q_old)q(8)=q_old                                                    
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(8),sign)
     if(q(8)<q_old)q(8)=q_old  

     lnods_cat(8,1,1)=ipoin
     lnods_cat(8,1,2)=lext2(5)
     lnods_cat(8,1,3)=lext2(6)
     lnods_cat(8,1,4)=lext2(7)
     lnods_cat(8,2,1)=ipoin
     lnods_cat(8,2,2)=lext2(1)
     lnods_cat(8,2,3)=lext2(5)
     lnods_cat(8,2,4)=lext2(7)
     lnods_cat(8,3,1)=ipoin
     lnods_cat(8,3,2)=lext2(1)
     lnods_cat(8,3,3)=lext2(3)
     lnods_cat(8,3,4)=lext2(5)
     lnods_cat(8,4,1)=ipoin
     lnods_cat(8,4,2)=lext2(3)
     lnods_cat(8,4,3)=lext2(4)
     lnods_cat(8,4,4)=lext2(5)
     lnods_cat(8,5,1)=ipoin
     lnods_cat(8,5,2)=lext2(1)
     lnods_cat(8,5,3)=lext2(2)
     lnods_cat(8,5,4)=lext2(3)


     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(9),sign)   
     q_old = q(9)                                                                
     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(6),q(9),sign)
     if(q(9)<q_old)q(9)=q_old                                                    
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(5),q(9),sign)
     if(q(9)<q_old)q(9)=q_old                                                    
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(9),sign)
     if(q(9)<q_old)q(9)=q_old                                                    
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(9),sign)
     if(q(9)<q_old)q(9)=q_old  

     lnods_cat(9,1,1)=ipoin
     lnods_cat(9,1,2)=lext2(1)
     lnods_cat(9,1,3)=lext2(6)
     lnods_cat(9,1,4)=lext2(7)
     lnods_cat(9,2,1)=ipoin
     lnods_cat(9,2,2)=lext2(1)
     lnods_cat(9,2,3)=lext2(5)
     lnods_cat(9,2,4)=lext2(6)
     lnods_cat(9,3,1)=ipoin
     lnods_cat(9,3,2)=lext2(1)
     lnods_cat(9,3,3)=lext2(3)
     lnods_cat(9,3,4)=lext2(5)
     lnods_cat(9,4,1)=ipoin
     lnods_cat(9,4,2)=lext2(3)
     lnods_cat(9,4,3)=lext2(4)
     lnods_cat(9,4,4)=lext2(5)
     lnods_cat(9,5,1)=ipoin
     lnods_cat(9,5,2)=lext2(1)
     lnods_cat(9,5,3)=lext2(2)
     lnods_cat(9,5,4)=lext2(3)

     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(10),sign)   
     q_old = q(10)                                                                
     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(7),q(10),sign)
     if(q(10)<q_old)q(10)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(5),q(10),sign)
     if(q(10)<q_old)q(10)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(5),q(10),sign)
     if(q(10)<q_old)q(10)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(10),sign)
     if(q(10)<q_old)q(10)=q_old  

     lnods_cat(10,1,1)=ipoin
     lnods_cat(10,1,2)=lext2(5)
     lnods_cat(10,1,3)=lext2(6)
     lnods_cat(10,1,4)=lext2(7)
     lnods_cat(10,2,1)=ipoin
     lnods_cat(10,2,2)=lext2(1)
     lnods_cat(10,2,3)=lext2(5)
     lnods_cat(10,2,4)=lext2(7)
     lnods_cat(10,3,1)=ipoin
     lnods_cat(10,3,2)=lext2(1)
     lnods_cat(10,3,3)=lext2(2)
     lnods_cat(10,3,4)=lext2(5)
     lnods_cat(10,4,1)=ipoin
     lnods_cat(10,4,2)=lext2(2)
     lnods_cat(10,4,3)=lext2(3)
     lnods_cat(10,4,4)=lext2(5)
     lnods_cat(10,5,1)=ipoin
     lnods_cat(10,5,2)=lext2(3)
     lnods_cat(10,5,3)=lext2(4)
     lnods_cat(10,5,4)=lext2(5)

     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(11),sign)   
     q_old = q(11)                                                                
     call dod_qualit(ipoin,lext2(2),lext2(5),lext2(7),q(11),sign)
     if(q(11)<q_old)q(11)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(11),sign)
     if(q(11)<q_old)q(11)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(5),q(11),sign)
     if(q(11)<q_old)q(11)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(11),sign)
     if(q(11)<q_old)q(11)=q_old  
     lnods_cat(11,1,1)=ipoin
     lnods_cat(11,1,2)=lext2(5)
     lnods_cat(11,1,3)=lext2(6)
     lnods_cat(11,1,4)=lext2(7)
     lnods_cat(11,2,1)=ipoin
     lnods_cat(11,2,2)=lext2(2)
     lnods_cat(11,2,3)=lext2(5)
     lnods_cat(11,2,4)=lext2(7)
     lnods_cat(11,3,1)=ipoin
     lnods_cat(11,3,2)=lext2(1)
     lnods_cat(11,3,3)=lext2(2)
     lnods_cat(11,3,4)=lext2(7)
     lnods_cat(11,4,1)=ipoin
     lnods_cat(11,4,2)=lext2(2)
     lnods_cat(11,4,3)=lext2(3)
     lnods_cat(11,4,4)=lext2(5)
     lnods_cat(11,5,1)=ipoin
     lnods_cat(11,5,2)=lext2(3)
     lnods_cat(11,5,3)=lext2(4)
     lnods_cat(11,5,4)=lext2(5)


     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(12),sign)   
     q_old = q(12)                                                                
     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(6),q(12),sign)
     if(q(12)<q_old)q(12)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(5),q(12),sign)
     if(q(12)<q_old)q(12)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(5),q(12),sign)
     if(q(12)<q_old)q(12)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(12),sign)
     if(q(12)<q_old)q(12)=q_old  
     lnods_cat(12,1,1)=ipoin
     lnods_cat(12,1,2)=lext2(1)
     lnods_cat(12,1,3)=lext2(6)
     lnods_cat(12,1,4)=lext2(7)
     lnods_cat(12,2,1)=ipoin
     lnods_cat(12,2,2)=lext2(1)
     lnods_cat(12,2,3)=lext2(5)
     lnods_cat(12,2,4)=lext2(6)
     lnods_cat(12,3,1)=ipoin
     lnods_cat(12,3,2)=lext2(1)
     lnods_cat(12,3,3)=lext2(2)
     lnods_cat(12,3,4)=lext2(5)
     lnods_cat(12,4,1)=ipoin
     lnods_cat(12,4,2)=lext2(2)
     lnods_cat(12,4,3)=lext2(3)
     lnods_cat(12,4,4)=lext2(5)
     lnods_cat(12,5,1)=ipoin
     lnods_cat(12,5,2)=lext2(3)
     lnods_cat(12,5,3)=lext2(4)
     lnods_cat(12,5,4)=lext2(5)


     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(13),sign)   
     q_old = q(13)                                                                
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(6),q(13),sign)
     if(q(13)<q_old)q(13)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(5),lext2(6),q(13),sign)
     if(q(13)<q_old)q(13)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(5),q(13),sign)
     if(q(13)<q_old)q(13)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(13),sign)
     if(q(13)<q_old)q(13)=q_old  
     lnods_cat(13,1,1)=ipoin
     lnods_cat(13,1,2)=lext2(5)
     lnods_cat(13,1,3)=lext2(6)
     lnods_cat(13,1,4)=lext2(7)
     lnods_cat(13,2,1)=ipoin
     lnods_cat(13,2,2)=lext2(1)
     lnods_cat(13,2,3)=lext2(2)
     lnods_cat(13,2,4)=lext2(6)
     lnods_cat(13,3,1)=ipoin
     lnods_cat(13,3,2)=lext2(2)
     lnods_cat(13,3,3)=lext2(5)
     lnods_cat(13,3,4)=lext2(6)
     lnods_cat(13,4,1)=ipoin
     lnods_cat(13,4,2)=lext2(2)
     lnods_cat(13,4,3)=lext2(3)
     lnods_cat(13,4,4)=lext2(5)
     lnods_cat(13,5,1)=ipoin
     lnods_cat(13,5,2)=lext2(3)
     lnods_cat(13,5,3)=lext2(4)
     lnods_cat(13,5,4)=lext2(5)

     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(14),sign)   
     q_old = q(14)                                                                
     call dod_qualit(ipoin,lext2(2),lext2(6),lext2(7),q(14),sign)
     if(q(14)<q_old)q(14)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(5),lext2(6),q(14),sign)
     if(q(14)<q_old)q(14)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(5),q(14),sign)
     if(q(14)<q_old)q(14)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(5),q(14),sign)
     if(q(14)<q_old)q(14)=q_old  
     lnods_cat(14,1,1)=ipoin
     lnods_cat(14,1,2)=lext2(1)
     lnods_cat(14,1,3)=lext2(2)
     lnods_cat(14,1,4)=lext2(7)
     lnods_cat(14,2,1)=ipoin
     lnods_cat(14,2,2)=lext2(2)
     lnods_cat(14,2,3)=lext2(6)
     lnods_cat(14,2,4)=lext2(7)
     lnods_cat(14,3,1)=ipoin
     lnods_cat(14,3,2)=lext2(2)
     lnods_cat(14,3,3)=lext2(5)
     lnods_cat(14,3,4)=lext2(6)
     lnods_cat(14,4,1)=ipoin
     lnods_cat(14,4,2)=lext2(2)
     lnods_cat(14,4,3)=lext2(3)
     lnods_cat(14,4,4)=lext2(5)
     lnods_cat(14,5,1)=ipoin
     lnods_cat(14,5,2)=lext2(3)
     lnods_cat(14,5,3)=lext2(4)
     lnods_cat(14,5,4)=lext2(5)

     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(15),sign)   
     q_old = q(15)                                                                
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(7),q(15),sign)
     if(q(15)<q_old)q(15)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(6),lext2(7),q(15),sign)
     if(q(15)<q_old)q(15)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(6),q(15),sign)
     if(q(15)<q_old)q(15)=q_old                                                   
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(15),sign)
     if(q(15)<q_old)q(15)=q_old  
     lnods_cat(15,1,1)=ipoin
     lnods_cat(15,1,2)=lext2(1)
     lnods_cat(15,1,3)=lext2(2)
     lnods_cat(15,1,4)=lext2(3)
     lnods_cat(15,2,1)=ipoin
     lnods_cat(15,2,2)=lext2(1)
     lnods_cat(15,2,3)=lext2(3)
     lnods_cat(15,2,4)=lext2(7)
     lnods_cat(15,3,1)=ipoin
     lnods_cat(15,3,2)=lext2(3)
     lnods_cat(15,3,3)=lext2(6)
     lnods_cat(15,3,4)=lext2(7)
     lnods_cat(15,4,1)=ipoin
     lnods_cat(15,4,2)=lext2(3)
     lnods_cat(15,4,3)=lext2(4)
     lnods_cat(15,4,4)=lext2(6)
     lnods_cat(15,5,1)=ipoin
     lnods_cat(15,5,2)=lext2(4)
     lnods_cat(15,5,3)=lext2(5)
     lnods_cat(15,5,4)=lext2(6)

     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(16),sign)   
     q_old = q(16)                                                                
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(7),q(16),sign)
     if(q(16)<q_old)q(16)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(7),q(16),sign)
     if(q(16)<q_old)q(16)=q_old                                                   
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(7),q(16),sign)
     if(q(16)<q_old)q(16)=q_old                                                   
     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(16),sign)
     if(q(16)<q_old)q(16)=q_old  
     lnods_cat(16,1,1)=ipoin
     lnods_cat(16,1,2)=lext2(1)
     lnods_cat(16,1,3)=lext2(2)
     lnods_cat(16,1,4)=lext2(3)
     lnods_cat(16,2,1)=ipoin
     lnods_cat(16,2,2)=lext2(2)
     lnods_cat(16,2,3)=lext2(3)
     lnods_cat(16,2,4)=lext2(7)
     lnods_cat(16,3,1)=ipoin
     lnods_cat(16,3,2)=lext2(3)
     lnods_cat(16,3,3)=lext2(4)
     lnods_cat(16,3,4)=lext2(7)
     lnods_cat(16,4,1)=ipoin
     lnods_cat(16,4,2)=lext2(4)
     lnods_cat(16,4,3)=lext2(5)
     lnods_cat(16,4,4)=lext2(7)
     lnods_cat(16,5,1)=ipoin
     lnods_cat(16,5,2)=lext2(5)
     lnods_cat(16,5,3)=lext2(6)
     lnods_cat(16,5,4)=lext2(7)

     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(17),sign)   
     q_old = q(17)                                                                
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(6),q(17),sign)
     if(q(17)<q_old)q(17)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(6),q(17),sign)
     if(q(17)<q_old)q(17)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(17),sign)
     if(q(17)<q_old)q(17)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(17),sign)
     if(q(17)<q_old)q(17)=q_old  
     lnods_cat(17,1,1)=ipoin
     lnods_cat(17,1,2)=lext2(4)
     lnods_cat(17,1,3)=lext2(5)
     lnods_cat(17,1,4)=lext2(6)
     lnods_cat(17,2,1)=ipoin
     lnods_cat(17,2,2)=lext2(3)
     lnods_cat(17,2,3)=lext2(4)
     lnods_cat(17,2,4)=lext2(6)
     lnods_cat(17,3,1)=ipoin
     lnods_cat(17,3,2)=lext2(1)
     lnods_cat(17,3,3)=lext2(3)
     lnods_cat(17,3,4)=lext2(6)
     lnods_cat(17,4,1)=ipoin
     lnods_cat(17,4,2)=lext2(1)
     lnods_cat(17,4,3)=lext2(2)
     lnods_cat(17,4,4)=lext2(3)
     lnods_cat(17,5,1)=ipoin
     lnods_cat(17,5,2)=lext2(1)
     lnods_cat(17,5,3)=lext2(6)
     lnods_cat(17,5,4)=lext2(7)


     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(18),sign)   
     q_old = q(18)                                                                
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(6),q(18),sign)
     if(q(18)<q_old)q(18)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(6),q(18),sign)
     if(q(18)<q_old)q(18)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(6),q(18),sign)
     if(q(18)<q_old)q(18)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(18),sign)
     if(q(18)<q_old)q(18)=q_old  

     lnods_cat(18,1,1)=ipoin
     lnods_cat(18,1,2)=lext2(4)
     lnods_cat(18,1,3)=lext2(5)
     lnods_cat(18,1,4)=lext2(6)
     lnods_cat(18,2,1)=ipoin
     lnods_cat(18,2,2)=lext2(3)
     lnods_cat(18,2,3)=lext2(4)
     lnods_cat(18,2,4)=lext2(6)
     lnods_cat(18,3,1)=ipoin
     lnods_cat(18,3,2)=lext2(2)
     lnods_cat(18,3,3)=lext2(3)
     lnods_cat(18,3,4)=lext2(6)
     lnods_cat(18,4,1)=ipoin
     lnods_cat(18,4,2)=lext2(1)
     lnods_cat(18,4,3)=lext2(2)
     lnods_cat(18,4,4)=lext2(6)
     lnods_cat(18,5,1)=ipoin
     lnods_cat(18,5,2)=lext2(1)
     lnods_cat(18,5,3)=lext2(6)
     lnods_cat(18,5,4)=lext2(7)

     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(19),sign)   
     q_old = q(19)                                                                
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(6),q(19),sign)
     if(q(19)<q_old)q(19)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(6),q(19),sign)
     if(q(19)<q_old)q(19)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(6),lext2(7),q(19),sign)
     if(q(19)<q_old)q(19)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(19),sign)
     if(q(19)<q_old)q(19)=q_old  

     lnods_cat(19,1,1)=ipoin
     lnods_cat(19,1,2)=lext2(4)
     lnods_cat(19,1,3)=lext2(5)
     lnods_cat(19,1,4)=lext2(6)
     lnods_cat(19,2,1)=ipoin
     lnods_cat(19,2,2)=lext2(3)
     lnods_cat(19,2,3)=lext2(4)
     lnods_cat(19,2,4)=lext2(6)
     lnods_cat(19,3,1)=ipoin
     lnods_cat(19,3,2)=lext2(2)
     lnods_cat(19,3,3)=lext2(3)
     lnods_cat(19,3,4)=lext2(6)
     lnods_cat(19,4,1)=ipoin
     lnods_cat(19,4,2)=lext2(2)
     lnods_cat(19,4,3)=lext2(6)
     lnods_cat(19,4,4)=lext2(7)
     lnods_cat(19,5,1)=ipoin
     lnods_cat(19,5,2)=lext2(1)
     lnods_cat(19,5,3)=lext2(2)
     lnods_cat(19,5,4)=lext2(7)

     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(20),sign)   
     q_old = q(20)                                                                
     call dod_qualit(ipoin,lext2(4),lext2(6),lext2(7),q(20),sign)
     if(q(20)<q_old)q(20)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(7),q(20),sign)
     if(q(20)<q_old)q(20)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(7),q(20),sign)
     if(q(20)<q_old)q(20)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(20),sign)
     if(q(20)<q_old)q(20)=q_old  

     lnods_cat(20,1,1)=ipoin
     lnods_cat(20,1,2)=lext2(4)
     lnods_cat(20,1,3)=lext2(5)
     lnods_cat(20,1,4)=lext2(6)
     lnods_cat(20,2,1)=ipoin
     lnods_cat(20,2,2)=lext2(4)
     lnods_cat(20,2,3)=lext2(6)
     lnods_cat(20,2,4)=lext2(7)
     lnods_cat(20,3,1)=ipoin
     lnods_cat(20,3,2)=lext2(3)
     lnods_cat(20,3,3)=lext2(4)
     lnods_cat(20,3,4)=lext2(7)
     lnods_cat(20,4,1)=ipoin
     lnods_cat(20,4,2)=lext2(1)
     lnods_cat(20,4,3)=lext2(3)
     lnods_cat(20,4,4)=lext2(7)
     lnods_cat(20,5,1)=ipoin
     lnods_cat(20,5,2)=lext2(1)
     lnods_cat(20,5,3)=lext2(2)
     lnods_cat(20,5,4)=lext2(3)

     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(21),sign)   
     q_old = q(21)                                                                
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(7),q(21),sign)
     if(q(21)<q_old)q(21)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(7),q(21),sign)
     if(q(21)<q_old)q(1)=q_old                                                    
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(7),q(21),sign)
     if(q(21)<q_old)q(21)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(21),sign)
     if(q(21)<q_old)q(21)=q_old  
     lnods_cat(21,1,1)=ipoin
     lnods_cat(21,1,2)=lext2(5)
     lnods_cat(21,1,3)=lext2(6)
     lnods_cat(21,1,4)=lext2(7)
     lnods_cat(21,2,1)=ipoin
     lnods_cat(21,2,2)=lext2(4)
     lnods_cat(21,2,3)=lext2(5)
     lnods_cat(21,2,4)=lext2(7)
     lnods_cat(21,3,1)=ipoin
     lnods_cat(21,3,2)=lext2(3)
     lnods_cat(21,3,3)=lext2(4)
     lnods_cat(21,3,4)=lext2(7)
     lnods_cat(21,4,1)=ipoin
     lnods_cat(21,4,2)=lext2(1)
     lnods_cat(21,4,3)=lext2(3)
     lnods_cat(21,4,4)=lext2(7)
     lnods_cat(21,5,1)=ipoin
     lnods_cat(21,5,2)=lext2(1)
     lnods_cat(21,5,3)=lext2(2)
     lnods_cat(21,5,4)=lext2(3)


     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(22),sign)   
     q_old = q(22)                                                                
     call dod_qualit(ipoin,lext2(4),lext2(6),lext2(7),q(22),sign)
     if(q(22)<q_old)q(22)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(7),q(22),sign)
     if(q(22)<q_old)q(22)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(7),q(22),sign)
     if(q(22)<q_old)q(22)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(22),sign)
     if(q(22)<q_old)q(22)=q_old  
     lnods_cat(22,1,1)=ipoin
     lnods_cat(22,1,2)=lext2(4)
     lnods_cat(22,1,3)=lext2(5)
     lnods_cat(22,1,4)=lext2(6)
     lnods_cat(22,2,1)=ipoin
     lnods_cat(22,2,2)=lext2(4)
     lnods_cat(22,2,3)=lext2(6)
     lnods_cat(22,2,4)=lext2(7)
     lnods_cat(22,3,1)=ipoin
     lnods_cat(22,3,2)=lext2(3)
     lnods_cat(22,3,3)=lext2(4)
     lnods_cat(22,3,4)=lext2(7)
     lnods_cat(22,4,1)=ipoin
     lnods_cat(22,4,2)=lext2(2)
     lnods_cat(22,4,3)=lext2(3)
     lnods_cat(22,4,4)=lext2(7)
     lnods_cat(22,5,1)=ipoin
     lnods_cat(22,5,2)=lext2(1)
     lnods_cat(22,5,3)=lext2(2)
     lnods_cat(22,5,4)=lext2(7)

     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(23),sign)   
     q_old = q(23)                                                                
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(7),q(23),sign)
     if(q(23)<q_old)q(23)=q_old                                                   
     call dod_qualit(ipoin,lext2(3),lext2(4),lext2(7),q(23),sign)
     if(q(23)<q_old)q(23)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(7),q(23),sign)
     if(q(23)<q_old)q(23)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(23),sign)
     if(q(23)<q_old)q(23)=q_old  
     lnods_cat(23,1,1)=ipoin
     lnods_cat(23,1,2)=lext2(5)
     lnods_cat(23,1,3)=lext2(6)
     lnods_cat(23,1,4)=lext2(7)
     lnods_cat(23,2,1)=ipoin
     lnods_cat(23,2,2)=lext2(4)
     lnods_cat(23,2,3)=lext2(5)
     lnods_cat(23,2,4)=lext2(7)
     lnods_cat(23,3,1)=ipoin
     lnods_cat(23,3,2)=lext2(3)
     lnods_cat(23,3,3)=lext2(4)
     lnods_cat(23,3,4)=lext2(7)
     lnods_cat(23,4,1)=ipoin
     lnods_cat(23,4,2)=lext2(2)
     lnods_cat(23,4,3)=lext2(3)
     lnods_cat(23,4,4)=lext2(7)
     lnods_cat(23,5,1)=ipoin
     lnods_cat(23,5,2)=lext2(1)
     lnods_cat(23,5,3)=lext2(2)
     lnods_cat(23,5,4)=lext2(7)

     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(24),sign)   
     q_old = q(24)                                                                
     call dod_qualit(ipoin,lext2(4),lext2(6),lext2(7),q(24),sign)
     if(q(24)<q_old)q(24)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(7),q(24),sign)
     if(q(24)<q_old)q(24)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(4),q(24),sign)
     if(q(24)<q_old)q(24)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(24),sign)
     if(q(24)<q_old)q(24)=q_old  
     lnods_cat(24,1,1)=ipoin
     lnods_cat(24,1,2)=lext2(4)
     lnods_cat(24,1,3)=lext2(5)
     lnods_cat(24,1,4)=lext2(6)
     lnods_cat(24,2,1)=ipoin
     lnods_cat(24,2,2)=lext2(4)
     lnods_cat(24,2,3)=lext2(6)
     lnods_cat(24,2,4)=lext2(7)
     lnods_cat(24,3,1)=ipoin
     lnods_cat(24,3,2)=lext2(1)
     lnods_cat(24,3,3)=lext2(4)
     lnods_cat(24,3,4)=lext2(7)
     lnods_cat(24,4,1)=ipoin
     lnods_cat(24,4,2)=lext2(1)
     lnods_cat(24,4,3)=lext2(3)
     lnods_cat(24,4,4)=lext2(4)
     lnods_cat(24,5,1)=ipoin
     lnods_cat(24,5,2)=lext2(1)
     lnods_cat(24,5,3)=lext2(2)
     lnods_cat(24,5,4)=lext2(3)

     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(25),sign)   
     q_old = q(25)                                                                
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(6),q(25),sign)
     if(q(25)<q_old)q(25)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(4),q(25),sign)
     if(q(25)<q_old)q(25)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(25),sign)
     if(q(25)<q_old)q(25)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(25),sign)
     if(q(25)<q_old)q(25)=q_old  
     lnods_cat(25,1,1)=ipoin
     lnods_cat(25,1,2)=lext2(4)
     lnods_cat(25,1,3)=lext2(5)
     lnods_cat(25,1,4)=lext2(6)
     lnods_cat(25,2,1)=ipoin
     lnods_cat(25,2,2)=lext2(1)
     lnods_cat(25,2,3)=lext2(4)
     lnods_cat(25,2,4)=lext2(6)
     lnods_cat(25,3,1)=ipoin
     lnods_cat(25,3,2)=lext2(1)
     lnods_cat(25,3,3)=lext2(3)
     lnods_cat(25,3,4)=lext2(4)
     lnods_cat(25,4,1)=ipoin
     lnods_cat(25,4,2)=lext2(1)
     lnods_cat(25,4,3)=lext2(2)
     lnods_cat(25,4,4)=lext2(3)
     lnods_cat(25,5,1)=ipoin
     lnods_cat(25,5,2)=lext2(1)
     lnods_cat(25,5,3)=lext2(6)
     lnods_cat(25,5,4)=lext2(7)

     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(26),sign)   
     q_old = q(26)                                                                
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(7),q(26),sign)
     if(q(26)<q_old)q(26)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(7),q(26),sign)
     if(q(26)<q_old)q(26)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(4),q(26),sign)
     if(q(26)<q_old)q(26)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(26),sign)
     if(q(26)<q_old)q(26)=q_old  
     lnods_cat(26,1,1)=ipoin
     lnods_cat(26,1,2)=lext2(5)
     lnods_cat(26,1,3)=lext2(6)
     lnods_cat(26,1,4)=lext2(7)
     lnods_cat(26,2,1)=ipoin
     lnods_cat(26,2,2)=lext2(4)
     lnods_cat(26,2,3)=lext2(5)
     lnods_cat(26,2,4)=lext2(7)
     lnods_cat(26,3,1)=ipoin
     lnods_cat(26,3,2)=lext2(1)
     lnods_cat(26,3,3)=lext2(4)
     lnods_cat(26,3,4)=lext2(7)
     lnods_cat(26,4,1)=ipoin
     lnods_cat(26,4,2)=lext2(1)
     lnods_cat(26,4,3)=lext2(3)
     lnods_cat(26,4,4)=lext2(4)
     lnods_cat(26,5,1)=ipoin
     lnods_cat(26,5,2)=lext2(1)
     lnods_cat(26,5,3)=lext2(2)
     lnods_cat(26,5,4)=lext2(3)

     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(27),sign)
     if(q(27)<q_old)q(27)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(7),q(27),sign)
     if(q(27)<q_old)q(27)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(5),q(27),sign)
     if(q(27)<q_old)q(27)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(4),q(27),sign)
     if(q(27)<q_old)q(27)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(27),sign)
     if(q(27)<q_old)q(27)=q_old  
     lnods_cat(27,1,1)=ipoin
     lnods_cat(27,1,2)=lext2(5)
     lnods_cat(27,1,3)=lext2(6)
     lnods_cat(27,1,4)=lext2(7)
     lnods_cat(27,2,1)=ipoin
     lnods_cat(27,2,2)=lext2(1)
     lnods_cat(27,2,3)=lext2(5)
     lnods_cat(27,2,4)=lext2(7)
     lnods_cat(27,3,1)=ipoin
     lnods_cat(27,3,2)=lext2(1)
     lnods_cat(27,3,3)=lext2(4)
     lnods_cat(27,3,4)=lext2(5)
     lnods_cat(27,4,1)=ipoin
     lnods_cat(27,4,2)=lext2(1)
     lnods_cat(27,4,3)=lext2(3)
     lnods_cat(27,4,4)=lext2(4)
     lnods_cat(27,5,1)=ipoin
     lnods_cat(27,5,2)=lext2(1)
     lnods_cat(27,5,3)=lext2(2)
     lnods_cat(27,5,4)=lext2(3)

     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(28),sign)
     if(q(28)<q_old)q(28)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(6),q(28),sign)
     if(q(28)<q_old)q(28)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(5),q(28),sign)
     if(q(28)<q_old)q(28)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(3),lext2(4),q(28),sign)
     if(q(28)<q_old)q(28)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(3),q(28),sign)
     if(q(28)<q_old)q(28)=q_old  
     lnods_cat(28,1,1)=ipoin
     lnods_cat(28,1,2)=lext2(1)
     lnods_cat(28,1,3)=lext2(6)
     lnods_cat(28,1,4)=lext2(7)
     lnods_cat(28,2,1)=ipoin
     lnods_cat(28,2,2)=lext2(1)
     lnods_cat(28,2,3)=lext2(5)
     lnods_cat(28,2,4)=lext2(6)
     lnods_cat(28,3,1)=ipoin
     lnods_cat(28,3,2)=lext2(1)
     lnods_cat(28,3,3)=lext2(4)
     lnods_cat(28,3,4)=lext2(5)
     lnods_cat(28,4,1)=ipoin
     lnods_cat(28,4,2)=lext2(1)
     lnods_cat(28,4,3)=lext2(3)
     lnods_cat(28,4,4)=lext2(4)
     lnods_cat(28,5,1)=ipoin
     lnods_cat(28,5,2)=lext2(1)
     lnods_cat(28,5,3)=lext2(2)
     lnods_cat(28,5,4)=lext2(3)

     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(29),sign)
     if(q(29)<q_old)q(29)=q_old                                                   
     call dod_qualit(ipoin,lext2(4),lext2(6),lext2(7),q(29),sign)
     if(q(29)<q_old)q(29)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(7),q(29),sign)
     if(q(29)<q_old)q(29)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(4),q(29),sign)
     if(q(29)<q_old)q(29)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(29),sign)
     if(q(29)<q_old)q(29)=q_old  
     lnods_cat(29,1,1)=ipoin
     lnods_cat(29,1,2)=lext2(4)
     lnods_cat(29,1,3)=lext2(5)
     lnods_cat(29,1,4)=lext2(6)
     lnods_cat(29,2,1)=ipoin
     lnods_cat(29,2,2)=lext2(4)
     lnods_cat(29,2,3)=lext2(6)
     lnods_cat(29,2,4)=lext2(7)
     lnods_cat(29,3,1)=ipoin
     lnods_cat(29,3,2)=lext2(1)
     lnods_cat(29,3,3)=lext2(4)
     lnods_cat(29,3,4)=lext2(7)
     lnods_cat(29,4,1)=ipoin
     lnods_cat(29,4,2)=lext2(1)
     lnods_cat(29,4,3)=lext2(2)
     lnods_cat(29,4,4)=lext2(4)
     lnods_cat(29,5,1)=ipoin
     lnods_cat(29,5,2)=lext2(2)
     lnods_cat(29,5,3)=lext2(3)
     lnods_cat(29,5,4)=lext2(4)

     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(30),sign)
     if(q(30)<q_old)q(30)=q_old                                                   
     call dod_qualit(ipoin,lext2(4),lext2(6),lext2(7),q(30),sign)
     if(q(30)<q_old)q(30)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(7),q(30),sign)
     if(q(30)<q_old)q(30)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(30),sign)
     if(q(30)<q_old)q(30)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(30),sign)
     if(q(30)<q_old)q(30)=q_old  
     lnods_cat(30,1,1)=ipoin
     lnods_cat(30,1,2)=lext2(4)
     lnods_cat(30,1,3)=lext2(5)
     lnods_cat(30,1,4)=lext2(6)
     lnods_cat(30,2,1)=ipoin
     lnods_cat(30,2,2)=lext2(4)
     lnods_cat(30,2,3)=lext2(6)
     lnods_cat(30,2,4)=lext2(7)
     lnods_cat(30,3,1)=ipoin
     lnods_cat(30,3,2)=lext2(2)
     lnods_cat(30,3,3)=lext2(4)
     lnods_cat(30,3,4)=lext2(7)
     lnods_cat(30,4,1)=ipoin
     lnods_cat(30,4,2)=lext2(1)
     lnods_cat(30,4,3)=lext2(2)
     lnods_cat(30,4,4)=lext2(7)
     lnods_cat(30,5,1)=ipoin
     lnods_cat(30,5,2)=lext2(2)
     lnods_cat(30,5,3)=lext2(3)
     lnods_cat(30,5,4)=lext2(4)

     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(31),sign)
     if(q(31)<q_old)q(31)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(6),q(31),sign)
     if(q(31)<q_old)q(31)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(31),sign)
     if(q(31)<q_old)q(31)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(4),q(31),sign)
     if(q(31)<q_old)q(31)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(31),sign)
     if(q(31)<q_old)q(31)=q_old  
     lnods_cat(31,1,1)=ipoin
     lnods_cat(31,1,2)=lext2(4)
     lnods_cat(31,1,3)=lext2(5)
     lnods_cat(31,1,4)=lext2(6)
     lnods_cat(31,2,1)=ipoin
     lnods_cat(31,2,2)=lext2(1)
     lnods_cat(31,2,3)=lext2(4)
     lnods_cat(31,2,4)=lext2(6)
     lnods_cat(31,3,1)=ipoin
     lnods_cat(31,3,2)=lext2(1)
     lnods_cat(31,3,3)=lext2(6)
     lnods_cat(31,3,4)=lext2(7)
     lnods_cat(31,4,1)=ipoin
     lnods_cat(31,4,2)=lext2(1)
     lnods_cat(31,4,3)=lext2(2)
     lnods_cat(31,4,4)=lext2(4)
     lnods_cat(31,5,1)=ipoin
     lnods_cat(31,5,2)=lext2(2)
     lnods_cat(31,5,3)=lext2(3)
     lnods_cat(31,5,4)=lext2(4)

     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(32),sign)
     if(q(32)<q_old)q(32)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(32),sign)
     if(q(32)<q_old)q(32)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(6),q(32),sign)
     if(q(32)<q_old)q(32)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(6),q(32),sign)
     if(q(32)<q_old)q(32)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(32),sign)
     if(q(32)<q_old)q(32)=q_old  
     lnods_cat(32,1,1)=ipoin
     lnods_cat(32,1,2)=lext2(4)
     lnods_cat(32,1,3)=lext2(5)
     lnods_cat(32,1,4)=lext2(6)
     lnods_cat(32,2,1)=ipoin
     lnods_cat(32,2,2)=lext2(1)
     lnods_cat(32,2,3)=lext2(6)
     lnods_cat(32,2,4)=lext2(7)
     lnods_cat(32,3,1)=ipoin
     lnods_cat(32,3,2)=lext2(1)
     lnods_cat(32,3,3)=lext2(2)
     lnods_cat(32,3,4)=lext2(6)
     lnods_cat(32,4,1)=ipoin
     lnods_cat(32,4,2)=lext2(2)
     lnods_cat(32,4,3)=lext2(4)
     lnods_cat(32,4,4)=lext2(6)
     lnods_cat(32,5,1)=ipoin
     lnods_cat(32,5,2)=lext2(2)
     lnods_cat(32,5,3)=lext2(3)
     lnods_cat(32,5,4)=lext2(4)

     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(6),q(33),sign)
     if(q(33)<q_old)q(33)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(6),q(33),sign)
     if(q(33)<q_old)q(33)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(6),lext2(7),q(33),sign)
     if(q(33)<q_old)q(33)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(33),sign)
     if(q(33)<q_old)q(33)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(33),sign)
     if(q(33)<q_old)q(33)=q_old  
     lnods_cat(33,1,1)=ipoin
     lnods_cat(33,1,2)=lext2(4)
     lnods_cat(33,1,3)=lext2(5)
     lnods_cat(33,1,4)=lext2(6)
     lnods_cat(33,2,1)=ipoin
     lnods_cat(33,2,2)=lext2(2)
     lnods_cat(33,2,3)=lext2(4)
     lnods_cat(33,2,4)=lext2(6)
     lnods_cat(33,3,1)=ipoin
     lnods_cat(33,3,2)=lext2(2)
     lnods_cat(33,3,3)=lext2(6)
     lnods_cat(33,3,4)=lext2(7)
     lnods_cat(33,4,1)=ipoin
     lnods_cat(33,4,2)=lext2(1)
     lnods_cat(33,4,3)=lext2(2)
     lnods_cat(33,4,4)=lext2(7)
     lnods_cat(33,5,1)=ipoin
     lnods_cat(33,5,2)=lext2(2)
     lnods_cat(33,5,3)=lext2(3)
     lnods_cat(33,5,4)=lext2(4)

     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(34),sign)
     if(q(34)<q_old)q(34)=q_old                                                   
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(7),q(34),sign)
     if(q(34)<q_old)q(34)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(7),q(34),sign)
     if(q(34)<q_old)q(34)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(4),q(34),sign)
     if(q(34)<q_old)q(34)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(34),sign)
     if(q(34)<q_old)q(34)=q_old  
     lnods_cat(34,1,1)=ipoin
     lnods_cat(34,1,2)=lext2(5)
     lnods_cat(34,1,3)=lext2(6)
     lnods_cat(34,1,4)=lext2(7)
     lnods_cat(34,2,1)=ipoin
     lnods_cat(34,2,2)=lext2(4)
     lnods_cat(34,2,3)=lext2(5)
     lnods_cat(34,2,4)=lext2(7)
     lnods_cat(34,3,1)=ipoin
     lnods_cat(34,3,2)=lext2(1)
     lnods_cat(34,3,3)=lext2(4)
     lnods_cat(34,3,4)=lext2(7)
     lnods_cat(34,4,1)=ipoin
     lnods_cat(34,4,2)=lext2(1)
     lnods_cat(34,4,3)=lext2(2)
     lnods_cat(34,4,4)=lext2(4)
     lnods_cat(34,5,1)=ipoin
     lnods_cat(34,5,2)=lext2(2)
     lnods_cat(34,5,3)=lext2(3)
     lnods_cat(34,5,4)=lext2(4)

     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(35),sign)
     if(q(35)<q_old)q(35)=q_old                                                   
     call dod_qualit(ipoin,lext2(4),lext2(5),lext2(7),q(35),sign)
     if(q(35)<q_old)q(35)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(7),q(35),sign)
     if(q(35)<q_old)q(35)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(35),sign)
     if(q(35)<q_old)q(35)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(35),sign)
     if(q(35)<q_old)q(35)=q_old  
     lnods_cat(35,1,1)=ipoin
     lnods_cat(35,1,2)=lext2(5)
     lnods_cat(35,1,3)=lext2(6)
     lnods_cat(35,1,4)=lext2(7)
     lnods_cat(35,2,1)=ipoin
     lnods_cat(35,2,2)=lext2(4)
     lnods_cat(35,2,3)=lext2(5)
     lnods_cat(35,2,4)=lext2(7)
     lnods_cat(35,3,1)=ipoin
     lnods_cat(35,3,2)=lext2(2)
     lnods_cat(35,3,3)=lext2(4)
     lnods_cat(35,3,4)=lext2(7)
     lnods_cat(35,4,1)=ipoin
     lnods_cat(35,4,2)=lext2(2)
     lnods_cat(35,4,3)=lext2(3)
     lnods_cat(35,4,4)=lext2(4)
     lnods_cat(35,5,1)=ipoin
     lnods_cat(35,5,2)=lext2(1)
     lnods_cat(35,5,3)=lext2(2)
     lnods_cat(35,5,4)=lext2(7)

     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(36),sign)
     if(q(36)<q_old)q(36)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(7),q(36),sign)
     if(q(36)<q_old)q(36)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(5),q(36),sign)
     if(q(36)<q_old)q(36)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(4),q(36),sign)
     if(q(36)<q_old)q(36)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(36),sign)
     if(q(36)<q_old)q(36)=q_old  
     lnods_cat(36,1,1)=ipoin
     lnods_cat(36,1,2)=lext2(5)
     lnods_cat(36,1,3)=lext2(6)
     lnods_cat(36,1,4)=lext2(7)
     lnods_cat(36,2,1)=ipoin
     lnods_cat(36,2,2)=lext2(1)
     lnods_cat(36,2,3)=lext2(5)
     lnods_cat(36,2,4)=lext2(7)
     lnods_cat(36,3,1)=ipoin
     lnods_cat(36,3,2)=lext2(1)
     lnods_cat(36,3,3)=lext2(4)
     lnods_cat(36,3,4)=lext2(5)
     lnods_cat(36,4,1)=ipoin
     lnods_cat(36,4,2)=lext2(1)
     lnods_cat(36,4,3)=lext2(2)
     lnods_cat(36,4,4)=lext2(4)
     lnods_cat(36,5,1)=ipoin
     lnods_cat(36,5,2)=lext2(2)
     lnods_cat(36,5,3)=lext2(3)
     lnods_cat(36,5,4)=lext2(4)

     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(37),sign)
     if(q(37)<q_old)q(37)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(6),q(37),sign)
     if(q(37)<q_old)q(37)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(4),lext2(5),q(37),sign)
     if(q(37)<q_old)q(37)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(4),q(37),sign)
     if(q(37)<q_old)q(37)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(37),sign)
     if(q(37)<q_old)q(37)=q_old  
     lnods_cat(37,1,1)=ipoin
     lnods_cat(37,1,2)=lext2(1)
     lnods_cat(37,1,3)=lext2(6)
     lnods_cat(37,1,4)=lext2(7)
     lnods_cat(37,2,1)=ipoin
     lnods_cat(37,2,2)=lext2(1)
     lnods_cat(37,2,3)=lext2(5)
     lnods_cat(37,2,4)=lext2(6)
     lnods_cat(37,3,1)=ipoin
     lnods_cat(37,3,2)=lext2(1)
     lnods_cat(37,3,3)=lext2(4)
     lnods_cat(37,3,4)=lext2(5)
     lnods_cat(37,4,1)=ipoin
     lnods_cat(37,4,2)=lext2(1)
     lnods_cat(37,4,3)=lext2(2)
     lnods_cat(37,4,4)=lext2(4)
     lnods_cat(37,5,1)=ipoin
     lnods_cat(37,5,2)=lext2(2)
     lnods_cat(37,5,3)=lext2(3)
     lnods_cat(37,5,4)=lext2(4)

     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(38),sign)
     if(q(38)<q_old)q(38)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(7),q(38),sign)
     if(q(38)<q_old)q(38)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(5),q(38),sign)
     if(q(38)<q_old)q(38)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(5),q(38),sign)
     if(q(38)<q_old)q(38)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(38),sign)
     if(q(38)<q_old)q(38)=q_old  
     lnods_cat(38,1,1)=ipoin
     lnods_cat(38,1,2)=lext2(5)
     lnods_cat(38,1,3)=lext2(6)
     lnods_cat(38,1,4)=lext2(7)
     lnods_cat(38,2,1)=ipoin
     lnods_cat(38,2,2)=lext2(1)
     lnods_cat(38,2,3)=lext2(5)
     lnods_cat(38,2,4)=lext2(7)
     lnods_cat(38,3,1)=ipoin
     lnods_cat(38,3,2)=lext2(1)
     lnods_cat(38,3,3)=lext2(2)
     lnods_cat(38,3,4)=lext2(5)
     lnods_cat(38,4,1)=ipoin
     lnods_cat(38,4,2)=lext2(2)
     lnods_cat(38,4,3)=lext2(4)
     lnods_cat(38,4,4)=lext2(5)
     lnods_cat(38,5,1)=ipoin
     lnods_cat(38,5,2)=lext2(2)
     lnods_cat(38,5,3)=lext2(3)
     lnods_cat(38,5,4)=lext2(4)


     call dod_qualit(ipoin,lext2(5),lext2(6),lext2(7),q(39),sign)
     if(q(39)<q_old)q(39)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(5),lext2(7),q(39),sign)
     if(q(39)<q_old)q(39)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(5),q(39),sign)
     if(q(39)<q_old)q(39)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(39),sign)
     if(q(39)<q_old)q(39)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(39),sign)
     if(q(39)<q_old)q(39)=q_old  
     lnods_cat(39,1,1)=ipoin
     lnods_cat(39,1,2)=lext2(5)
     lnods_cat(39,1,3)=lext2(6)
     lnods_cat(39,1,4)=lext2(7)
     lnods_cat(39,2,1)=ipoin
     lnods_cat(39,2,2)=lext2(2)
     lnods_cat(39,2,3)=lext2(5)
     lnods_cat(39,2,4)=lext2(7)
     lnods_cat(39,3,1)=ipoin
     lnods_cat(39,3,2)=lext2(2)
     lnods_cat(39,3,3)=lext2(4)
     lnods_cat(39,3,4)=lext2(5)
     lnods_cat(39,4,1)=ipoin
     lnods_cat(39,4,2)=lext2(2)
     lnods_cat(39,4,3)=lext2(3)
     lnods_cat(39,4,4)=lext2(4)
     lnods_cat(39,5,1)=ipoin
     lnods_cat(39,5,2)=lext2(1)
     lnods_cat(39,5,3)=lext2(2)
     lnods_cat(39,5,4)=lext2(7)

     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(40),sign)
     if(q(40)<q_old)q(40)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(5),lext2(6),q(40),sign)
     if(q(40)<q_old)q(40)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(5),q(40),sign)
     if(q(40)<q_old)q(40)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(5),q(40),sign)
     if(q(40)<q_old)q(40)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(40),sign)
     if(q(40)<q_old)q(40)=q_old  
     lnods_cat(40,1,1)=ipoin
     lnods_cat(40,1,2)=lext2(1)
     lnods_cat(40,1,3)=lext2(6)
     lnods_cat(40,1,4)=lext2(7)
     lnods_cat(40,2,1)=ipoin
     lnods_cat(40,2,2)=lext2(1)
     lnods_cat(40,2,3)=lext2(5)
     lnods_cat(40,2,4)=lext2(6)
     lnods_cat(40,3,1)=ipoin
     lnods_cat(40,3,2)=lext2(1)
     lnods_cat(40,3,3)=lext2(2)
     lnods_cat(40,3,4)=lext2(5)
     lnods_cat(40,4,1)=ipoin
     lnods_cat(40,4,2)=lext2(2)
     lnods_cat(40,4,3)=lext2(4)
     lnods_cat(40,4,4)=lext2(5)
     lnods_cat(40,5,1)=ipoin
     lnods_cat(40,5,2)=lext2(2)
     lnods_cat(40,5,3)=lext2(3)
     lnods_cat(40,5,4)=lext2(4)

     call dod_qualit(ipoin,lext2(1),lext2(6),lext2(7),q(41),sign)
     if(q(41)<q_old)q(41)=q_old                                                   
     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(6),q(41),sign)
     if(q(41)<q_old)q(41)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(5),lext2(6),q(41),sign)
     if(q(41)<q_old)q(41)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(5),q(41),sign)
     if(q(41)<q_old)q(41)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(41),sign)
     if(q(41)<q_old)q(41)=q_old  
     lnods_cat(41,1,1)=ipoin
     lnods_cat(41,1,2)=lext2(1)
     lnods_cat(41,1,3)=lext2(6)
     lnods_cat(41,1,4)=lext2(7)
     lnods_cat(41,2,1)=ipoin
     lnods_cat(41,2,2)=lext2(1)
     lnods_cat(41,2,3)=lext2(2)
     lnods_cat(41,2,4)=lext2(6)
     lnods_cat(41,3,1)=ipoin
     lnods_cat(41,3,2)=lext2(2)
     lnods_cat(41,3,3)=lext2(5)
     lnods_cat(41,3,4)=lext2(6)
     lnods_cat(41,4,1)=ipoin
     lnods_cat(41,4,2)=lext2(2)
     lnods_cat(41,4,3)=lext2(4)
     lnods_cat(41,4,4)=lext2(5)
     lnods_cat(41,5,1)=ipoin
     lnods_cat(41,5,2)=lext2(2)
     lnods_cat(41,5,3)=lext2(3)
     lnods_cat(41,5,4)=lext2(4)

     call dod_qualit(ipoin,lext2(1),lext2(2),lext2(7),q(42),sign)
     if(q(42)<q_old)q(42)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(6),lext2(7),q(42),sign)
     if(q(42)<q_old)q(42)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(5),lext2(6),q(42),sign)
     if(q(42)<q_old)q(42)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(4),lext2(5),q(42),sign)
     if(q(42)<q_old)q(42)=q_old                                                   
     call dod_qualit(ipoin,lext2(2),lext2(3),lext2(4),q(42),sign)
     if(q(42)<q_old)q(42)=q_old  
     lnods_cat(42,1,1)=ipoin
     lnods_cat(42,1,2)=lext2(1)
     lnods_cat(42,1,3)=lext2(2)
     lnods_cat(42,1,4)=lext2(7)
     lnods_cat(42,2,1)=ipoin
     lnods_cat(42,2,2)=lext2(2)
     lnods_cat(42,2,3)=lext2(6)
     lnods_cat(42,2,4)=lext2(7)
     lnods_cat(42,3,1)=ipoin
     lnods_cat(42,3,2)=lext2(2)
     lnods_cat(42,3,3)=lext2(5)
     lnods_cat(42,3,4)=lext2(6)
     lnods_cat(42,4,1)=ipoin
     lnods_cat(42,4,2)=lext2(2)
     lnods_cat(42,4,3)=lext2(4)
     lnods_cat(42,4,4)=lext2(5)
     lnods_cat(42,5,1)=ipoin
     lnods_cat(42,5,2)=lext2(2)
     lnods_cat(42,5,3)=lext2(3)
     lnods_cat(42,5,4)=lext2(4)

     index = 1
     minim = q(1)

     do ii = 2,nucat
        if(q(ii)<minim ) then
           index = ii
           minim = q(ii)
        end if
     end do
     do itri = 1, nutri
!!!ielem_aux = -lntyp_dod(ipoin)%l(2*tboun+itri)
        ielem_aux = ielem_aux+1
        poex1=lnods_cat(index,itri,1)
        poex2=lnods_cat(index,itri,2)
        poex3=lnods_cat(index,itri,3)
        poex4=lnods_cat(index,itri,4)
        imatn          = lmatn_dod(poex1)
        ipoiz          = lpoiz_dod(poex1)
        lpext(number_fringe_nodes) % ltype(nboun_local+ielem_aux) = 30
        lpext(number_fringe_nodes) % lmate(nboun_local+ielem_aux) = imatn
        lpext(number_fringe_nodes) % lelez(nboun_local+ielem_aux) = ipoiz
        lpext(number_fringe_nodes) % lesub(nboun_local+ielem_aux) = isubd
        call dod_extens_qual3d(3_ip,poex1,poex2,poex3,poex4,kappa,asrad,dummr,sign,dumm5)
        !call dod_quaux(poex1,poex2,poex3,poex4)

        lpext(number_fringe_nodes) % lnods(1,nboun_local+ielem_aux)= poex1
        lpext(number_fringe_nodes) % lnods(2,nboun_local+ielem_aux)= poex2
        lpext(number_fringe_nodes) % lnods(3,nboun_local+ielem_aux)= poex3
        lpext(number_fringe_nodes) % lnods(4,nboun_local+ielem_aux)= poex4
        nelem_dod = nelem_dod + 1

        !print*,'concat5',lpext(number_fringe_nodes) % lnods(ielem_aux,:)
        if(lpext(number_fringe_nodes) % lnods(1,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(2,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(3,nboun_local+ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(4,nboun_local+ielem_aux)==0)then
           print*,'STOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOPPPPPP',ielem_aux,ipoin,nucat
           stop
        end if
     end do


  else
     print*,'EL NUMERO CATALAN ES MAYOR QUE 42!'
     stop


  end if
  call memchk(two,istat,mem_servi(1:2,servi),'LNODS_CAT','dod_concat',lnods_cat)
  deallocate(lnods_cat,stat=istat)
  if(istat/=0) call memerr(two,'LNODS_CAT','dod_concat',0_ip) 

  call memchk(two,istat,mem_servi(1:2,servi),'Q','dod_concat',q)
  deallocate(q,stat=istat)
  if(istat/=0) call memerr(two,'Q','dod_concat',0_ip)

  call memchk(two,istat,mem_servi(1:2,servi),'LEXT2_TMP','dod_concat',lext2_tmp)
  deallocate(lext2_tmp,stat=istat)
  if(istat/=0) call memerr(two,'LEXT"_TMP','dod_concat',0_ip) 

end subroutine dod_concat
