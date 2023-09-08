subroutine ibm_dpawal(iimbo,iwaib,coori,coorj,norma,dista) 
  !-----------------------------------------------------------------------
  ! NAME
  !    ibm_dpapar
  ! DESCRIPTION
  !    This routines find the minumium distance beetween two particles
  !    INPUT
  !       iimbo: id of the first particle
  !       jimbo: id of the second particle
  !    OUTPUT
  !       coori: point in iimbo nearest to jimbo
  !       coori: point in jimbo nearest to iimbo
  !       coori: number of faces of rthe particle
  !       dista: distance between particles
  !       norma: contact normal (use to solve collisions)
  ! USED BY
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  imbou,cutim,zeror
  use def_domain, only     :  ndime,nnode,mnoib,mnodb
  use def_immbou, only     :  twall_ibm
  use mod_kdtree

  implicit none
  integer(ip), intent(in)  :: iimbo,iwaib
  real(rp),    intent(out) :: dista,coori(3),coorj(3),norma(3)
  integer(ip)              :: idime,ipoib,iboib,inoib,jnoib,jboun,kboun
  real(rp)                 :: xcoor(3),tedis,tenor(3),proje(3),dummi
  real(rp)                 :: proji(3),projj(3),cooi1(3),cooi2(3)
  real(rp)                 :: bocod(ndime,mnoib),vect1(3),vect2(3),chkdi

  integer(ip)              :: indst,curr,xpoib,ilist,jlist
  real(rp)                 :: itemp,dismi1,dismi2,facto

  !
  ! Initialization
  !
  dista = 1.0e10_rp
  
  do ipoib = 1,imbou(iimbo) % npoib       
     xcoor(1)     = imbou(iimbo) % cooib(1,ipoib)
     xcoor(2)     = imbou(iimbo) % cooib(2,ipoib) 
     xcoor(ndime) = imbou(iimbo) % cooib(ndime,ipoib)

     twall_ibm(iwaib) % ldist(1) = 0.0_rp
     do idime = 1,ndime
        itemp  = max(xcoor(idime)-twall_ibm(iwaib) % sabox(1,idime,1), twall_ibm(iwaib) % sabox(2,idime,1)-xcoor(idime))
        twall_ibm(iwaib) % ldist(1) = twall_ibm(iwaib) % ldist(1) + itemp*itemp
     end do
     
     indst = 1_ip              
     twall_ibm(iwaib) % stru2(1) = 1_ip
                 
     do while( indst > 0_ip )
        curr  = twall_ibm(iwaib) % stru2(indst)
        indst = indst - 1_ip
        !
        ! Minimum distance
        !               
        if ( twall_ibm(iwaib) % ldist(curr) < dista ) then        

           if( twall_ibm(iwaib) % blink(curr) < 0_ip ) then              
              xpoib = ipoib   
              indst = 0_ip              
           else
              dismi1 = 0.0_rp
              dismi2 = 0.0_rp
              do idime = 1,ndime
                 itemp  = max(0.0_rp,twall_ibm(iwaib) % sabox(1,idime,twall_ibm(iwaib) % blink(curr)) - xcoor(idime))  
                 itemp  = max(itemp, xcoor(idime) - twall_ibm(iwaib) % sabox(2,idime,twall_ibm(iwaib) % blink(curr)))  
                 dismi1 = dismi1 + itemp * itemp
                 
                 itemp  = max(0.0_rp,twall_ibm(iwaib) % sabox(1,idime,twall_ibm(iwaib) % blink(curr)+1_ip) - xcoor(idime))  
                 itemp  = max(itemp, xcoor(idime) - twall_ibm(iwaib) % sabox(2,idime,twall_ibm(iwaib) % blink(curr)+1_ip))  
                 dismi2 = dismi2 + itemp * itemp
              end do
              indst = indst + 1_ip
              if (dismi1 > dismi2) then               
                 twall_ibm(iwaib) % ldist(twall_ibm(iwaib) % blink(curr)+1)  = dismi2
                 twall_ibm(iwaib) % stru2(indst)   = twall_ibm(iwaib) % blink(curr) + 1_ip
              else
                 twall_ibm(iwaib) % ldist(twall_ibm(iwaib) % blink(curr))    = dismi1             
                 twall_ibm(iwaib) % stru2(indst)   = twall_ibm(iwaib) % blink(curr)
              end if
           end if
        end if
     end do
  end do

  xcoor(1)     = imbou(iimbo) % cooib(1,xpoib)
  xcoor(2)     = imbou(iimbo) % cooib(2,xpoib)
  xcoor(ndime) = imbou(iimbo) % cooib(ndime,xpoib)
  
  call dpopar(1_ip,xcoor,&
       twall_ibm(iwaib) % npoin,mnoib,twall_ibm(iwaib) % nboun,1.0e10_rp,&
       twall_ibm(iwaib) % ltypb,twall_ibm(iwaib) % lnodb,twall_ibm(iwaib) % coord,&
       dista,norma,coorj,kboun,twall_ibm(iwaib) % fabox,twall_ibm(iwaib) % sabox,&
       twall_ibm(iwaib) % blink,twall_ibm(iwaib) % stru2,twall_ibm(iwaib) % ldist,&
       twall_ibm(iwaib) % lnele)

  chkdi = dista
  coori(1)     = xcoor(1)
  coori(2)     = xcoor(2)
  coori(ndime) = xcoor(ndime)

  !
  ! Determine the minimum distance beetween the jimbo NODES and the iimbo particle.
  !
  do ipoib = 1,twall_ibm(iwaib) % npoin

        xcoor(1)     = twall_ibm(iwaib) %coord(1,ipoib)
        xcoor(2)     = twall_ibm(iwaib) %coord(2,ipoib)
        xcoor(ndime) = twall_ibm(iwaib) %coord(ndime,ipoib)

        call dpopar(1_ip,xcoor,&
             imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,chkdi,&
             imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
             tedis,tenor,proje,kboun,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
             imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
             imbou(iimbo) % lnele)

        if ( tedis < dista ) then
           dista = tedis 
           chkdi = tedis
           
           coori(1)     = proje(1)
           coori(2)     = proje(2)
           coori(ndime) = proje(ndime)
           
           coorj(1)     = xcoor(1)
           coorj(2)     = xcoor(2)
           coorj(ndime) = xcoor(ndime)

           norma(1)     = -tenor(1)
           norma(2)     = -tenor(2)
           norma(ndime) = -tenor(ndime)              
        end if
  end do
  !
  ! Determine the minimum distance beetween the iimbo NODES and the jimbo particle
  !
  do ipoib = 1,imbou(iimbo) % npoib    

     xcoor(1)     = imbou(iimbo) % cooib(1,ipoib) 
     xcoor(2)     = imbou(iimbo) % cooib(2,ipoib) 
     xcoor(ndime) = imbou(iimbo) % cooib(ndime,ipoib)

     call dpopar(1_ip,xcoor,&
          twall_ibm(iwaib) % npoin,mnodb,twall_ibm(iwaib) % nboun,chkdi,&
          twall_ibm(iwaib) % ltypb,twall_ibm(iwaib) % lnodb,twall_ibm(iwaib) % coord,&
          tedis,tenor,proje,kboun,twall_ibm(iwaib) % fabox,twall_ibm(iwaib) % sabox,&
          twall_ibm(iwaib) % blink,twall_ibm(iwaib) % stru2,twall_ibm(iwaib) % ldist,&
          twall_ibm(iwaib) % lnele)

     if ( tedis < dista ) then
        dista = tedis 
        chkdi = tedis

        coori(1)     = xcoor(1)
        coori(2)     = xcoor(2)
        coori(ndime) = xcoor(ndime)

        coorj(1)     = proje(1)         
        coorj(2)     = proje(2)         
        coorj(ndime) = proje(ndime) 

        norma(1)     = tenor(1)
        norma(2)     = tenor(2)
        norma(ndime) = tenor(ndime)
     end if
 end do
  !
  ! Determine the minimum distance beetween the iimbo EDGES and the jimbo particle 
  !
  do iboib = 1,twall_ibm(iwaib) % nboun
        do inoib = 1,mnodb
           ipoib = twall_ibm(iwaib) % lnodb(inoib,iboib)            
           
           bocod(1    ,inoib) = twall_ibm(iwaib) %coord(1    ,ipoib)
           bocod(2    ,inoib) = twall_ibm(iwaib) %coord(2    ,ipoib)
           bocod(ndime,inoib) = twall_ibm(iwaib) %coord(ndime,ipoib)
        end do
        
        do inoib = 1,mnoib
           jnoib = inoib + 1
           if (jnoib > mnoib) jnoib = 1        
           call ibm_dsepar(iimbo,bocod(1,inoib),bocod(1,jnoib),chkdi,tedis,proji,projj,cooi1,cooi2)
           ! A new shorest distance have been found
           if ( tedis < dista ) then    
              dista = tedis 
              chkdi = tedis
              
              coori(1)     = proji(1)
              coori(2)     = proji(2)
              coori(ndime) = proji(ndime)

              coorj(1)     = projj(1)
              coorj(2)     = projj(2)
              coorj(ndime) = projj(ndime)       

              vect1(3) = 0.0_rp
              vect2(3) = 0.0_rp
              do idime = 1,ndime
                 vect1(idime) = bocod(idime,jnoib) - bocod(idime,inoib)                 
                 vect2(idime) = cooi1(idime) - cooi2(idime)
              end do
              call vecpro(vect1,vect2,norma,3)
              call vecuni(3_ip,norma,dummi)              

              tenor(1)     = 0.0_rp
              tenor(2)     = 0.0_rp
              tenor(ndime) = 0.0_rp
              do ilist =    1,twall_ibm(iwaib)%lnele(  twall_ibm(iwaib)%lnodb(inoib,iboib)  )%nelem
                 do jlist = 1,twall_ibm(iwaib)%lnele(  twall_ibm(iwaib)%lnodb(jnoib,iboib)  )%nelem
                    jboun =   twall_ibm(iwaib)%lnele(  twall_ibm(iwaib)%lnodb(inoib,iboib)  )%eleme(ilist)
                    kboun =   twall_ibm(iwaib)%lnele(  twall_ibm(iwaib)%lnodb(jnoib,iboib)  )%eleme(jlist)
                    if (jboun == kboun) then
                       call extbou(1_ip,nnode(twall_ibm(iwaib) %ltypb(jboun)),twall_ibm(iwaib) %lnodb(1,jboun),twall_ibm(iwaib) %coord,vect1)
                       do idime = 1,ndime  
                          tenor(idime) = tenor(idime) + vect1(idime)
                       end do
                       call vecuni(ndime,tenor,dummi)             
                    end if
                 end  do
              end  do

              facto = 0.0_rp
              do idime = 1,ndime
                 facto = facto +norma(idime)*tenor(idime)
              end do
              if (facto < 0.0_rp) then
                 do idime = 1,ndime
                    norma(idime) =  -norma(idime)
                 end do
              end if              
           end if
        end do
     end do

  if ( dista < -0.05_rp) then
          
     if ((iwaib==16)) then
        print *,dista,iimbo,iwaib
        print *,"jopeta total en las paredes!!!!"
     !   pause
     else
        dista = abs(dista)
     end if
  end if

end subroutine ibm_dpawal



