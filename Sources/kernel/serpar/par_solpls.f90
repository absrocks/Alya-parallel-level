subroutine par_solpls(itask)
  !------------------------------------------------------------------------
  !****f* Parall/par_solpls
  ! NAME
  !    par_sopls
  ! DESCRIPTION
  !    This routine sends graphs of boundary matrices 
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_parall
  use mod_memchk
  use mod_parall, only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall, only : commd,PAR_COMM_MY_CODE4
  use mod_parall, only : PAR_INTEGER
  use mod_parall, only : par_memor
  implicit none  
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,kpoin,jpoin,mpoin,kzdom,nnz
  integer(ip)             :: indice1,domai
  integer(4)              :: istat
  integer(ip), pointer    :: idx(:),pos(:),lnbin_loc(:)


  select case(itask)

  case(1)
     !
     ! Local-global numbering
     !
     if(kfl_paral==0) then
        !
        ! Master: send LNBIN_LOC
        !
        allocate( lnbin_loc(npoin_total),stat=istat)
        call memchk( zero, istat, par_memor, 'lnbin_loc', 'par_sengeo', lnbin_loc )
        do ipoin = 1, npoin_total
           jpoin            = lninv_loc(ipoin)
           lnbin_loc(ipoin) = lnper_par(jpoin)
        end do
        indice1 = 1
        do domai=1,npart_par
           npari =  npoin_par(domai)
           parin => lnbin_loc(indice1:)
           strin =  'LNBIN_LOC'
           call par_sendin()
           indice1 = indice1 + npoin_par(domai)
        end do
        call memchk( two, istat, par_memor, 'LNBIN_LOC','par_sengeo', lnbin_loc)
        deallocate( lnbin_loc, stat=istat )
        if(istat/=0) call memerr( two, 'LNBIN_LOC', 'par_sengeo',0_ip)

     else if(kfl_paral>0) then
        !
        ! Slave: receive LNBIN
        !
        npari =  npoin
        parin => lnbin
        call par_receiv()

     end if

  case(2)
     !
     ! Compute ABB matrix for each subdomain
     !
     if(kfl_paral==0) then
        !
        ! Master
        !
        allocate(idx(gnb+1),stat=istat)
        call memchk(zero,istat,par_memor,'WVERT_PAR','par_partit',idx)

        do domai=1,npart_par
           !
           ! Loop over subdomains
           !
           idx(1)=1
           do mpoin=1,ginde_par(2,domai)-1
              idx(mpoin+1)=1
           end do
           !
           ! Loop over own boundary
           !
           nnz=0
           do mpoin=ginde_par(2,domai),ginde_par(2,domai+1)-1
              !
              ! MPOIN: i/b boundary numbering (starting from 1)
              ! KPOIN: original numbering
              ! JPOIN: connected node to KPOIN in original mumbering
              !
              idx(mpoin+1)=idx(mpoin)
              kpoin=lninv_par(mpoin+gni)
              do kzdom=r_dom(kpoin),r_dom(kpoin+1)-1
                 jpoin=c_dom(kzdom)
                 if(lnper_par(jpoin)>gni) nnz=nnz+1
              end do
              idx(mpoin+1)=nnz+1
           end do
           !
           ! Rest of other boundary
           !
           do mpoin=ginde_par(2,domai+1)-1,gnb
              idx(mpoin+1)=nnz+1
           end do
           !
           ! Allocate POS
           !
           allocate(pos(nnz),stat=istat)
           call memchk(zero,istat,par_memor,'POS','par_solpls',pos)
           !
           ! Fill in POS
           !
           nnz=0
           do mpoin=ginde_par(2,domai),ginde_par(2,domai+1)-1
              !
              ! MPOIN: i/b boundary numbering (starting from 1)
              ! KPOIN: original numbering
              ! JPOIN: connected node to KPOIN in original mumbering
              !
              kpoin=lninv_par(mpoin+gni)
              do kzdom=r_dom(kpoin),r_dom(kpoin+1)-1
                 jpoin=c_dom(kzdom)
                 if(lnper_par(jpoin)>gni) then
                    nnz=nnz+1
                    pos(nnz)=lnper_par(jpoin)-gni
                 end if
              end do
           end do
           !
           ! Write matrix
           !
           !allocate(mat(nnz),stat=istat)
           !mat=1.0
           !ii=30+domai
           !call pspltm(&
           !  gnb,gnb,1_ip,0_ip,pos,idx,mat,&
           !  trim(title)//': '//naser(servi),0_ip,18.0_rp,'cm',&
           !  0_ip,0_ip,2_ip,ii)
           !deallocate(mat,stat=istat)

           !
           ! Send IDX and POS
           !    
           kfl_desti_par =  domai
           npari         =  gnb+1
           parin         => idx
           call par_sendin()
           npari         =  nnz
           parin         => pos
           call par_sendin()
           !
           ! Deallocate POS
           !
           call memchk(two,istat,par_memor,'POS','par_solpls',pos)
           deallocate(pos,stat=istat)
           if(istat/=0) call memerr(two,'POS','par_solpls',0_ip)

        end do
        !
        ! Deallocate IDX
        !
        call memchk(two,istat,par_memor,'IDX','par_solpls',idx)
        deallocate(idx,stat=istat)
        if(istat/=0) call memerr(two,'IDX','par_solpls',0_ip)

     else    
        
     end if

  end select

end subroutine par_solpls
