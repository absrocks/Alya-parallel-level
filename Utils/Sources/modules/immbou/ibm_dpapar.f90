subroutine ibm_dpapar(iimbo,jimbo,coori,coorj,norma,dista)         
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
  use def_master, only     :  imbou,cutim,zeror,ittim
  use def_domain, only     :  ndime,nnode,mnoib
  use mod_kdtree

  implicit none
  integer(ip), intent(in)    :: iimbo,jimbo
  real(rp),    intent(out )  :: dista,coori(3),coorj(3),norma(3)
  integer(ip)                :: idime,ipoib,iboib,jboib,kboib,inoib,jnoib,kboun
  real(rp)                   :: tedis,tenor(3),proje(3)
  real(rp)                   :: xcoor(3),chkdi
  real(rp)                   :: proji(3),projj(3),cooi1(3),cooi2(3)
  real(rp)                   :: bocod(ndime,mnoib),vect1(3),vect2(3)

  integer(ip)                :: indst,curr,xpoib,iboun,ilist,jlist
  real(rp)                   :: itemp,dismi1,dismi2,dummi,facto


  !
  ! Determine the minimum distance beetween the jimbo NODES and the iimbo particle.
  !  
  dista = 1.0e10_rp
  do ipoib = 1,imbou(jimbo) % npoib       
     xcoor(1)     = imbou(jimbo) % cooib(1,ipoib) 
     xcoor(2)     = imbou(jimbo) % cooib(2,ipoib) 
     xcoor(ndime) = imbou(jimbo) % cooib(ndime,ipoib)

     imbou(iimbo) % ldist(1) = 0.0_rp
     do idime = 1,ndime
        itemp  = max(xcoor(idime)-imbou(iimbo) % sabox(1,idime,1), imbou(iimbo) % sabox(2,idime,1)-xcoor(idime))
        imbou(iimbo) % ldist(1) = imbou(iimbo) % ldist(1) + itemp*itemp
     end do
     
     indst = 1_ip              
     imbou(iimbo) % struc(1) = 1_ip
                 
     do while( indst > 0_ip )
        curr  = imbou(iimbo) % struc(indst)
        indst = indst - 1_ip
        !
        ! Minimum distance
        !             
        if ( imbou(iimbo) % ldist(curr) < dista ) then        
           if( imbou(iimbo) % blink(curr) < 0_ip ) then              
              xpoib = ipoib   
              indst = 0_ip              
           else
              dismi1 = 0.0_rp
              dismi2 = 0.0_rp
              do idime = 1,ndime
                 itemp  = max(0.0_rp,imbou(iimbo) % sabox(1,idime,imbou(iimbo) % blink(curr)) - xcoor(idime))  
                 itemp  = max(itemp, xcoor(idime) - imbou(iimbo) % sabox(2,idime,imbou(iimbo) % blink(curr)))  
                 dismi1 = dismi1 + itemp * itemp
                 
                 itemp  = max(0.0_rp,imbou(iimbo) % sabox(1,idime,imbou(iimbo) % blink(curr)+1_ip) - xcoor(idime))  
                 itemp  = max(itemp, xcoor(idime) - imbou(iimbo) % sabox(2,idime,imbou(iimbo) % blink(curr)+1_ip))  
                 dismi2 = dismi2 + itemp * itemp
              end do
              indst = indst + 1_ip
              if (dismi1 > dismi2) then               
                 imbou(iimbo) % ldist(imbou(iimbo) % blink(curr)+1)  = dismi2
                 imbou(iimbo) % struc(indst)   = imbou(iimbo) % blink(curr) + 1_ip
              else
                 imbou(iimbo) % ldist(imbou(iimbo) % blink(curr))    = dismi1             
                 imbou(iimbo) % struc(indst)   = imbou(iimbo) % blink(curr)
              end if
           end if
        end if
     end do
  end do

  xcoor(1)     = imbou(jimbo) % cooib(1,xpoib) 
  xcoor(2)     = imbou(jimbo) % cooib(2,xpoib) 
  xcoor(ndime) = imbou(jimbo) % cooib(ndime,xpoib)
   
  call dpopar(1_ip,xcoor,&
       imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
       imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
       dista,norma,coori,kboun, imbou(iimbo) % fabox, imbou(iimbo) % sabox,&
       imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
       imbou(iimbo) % lnele)

  norma(1)     = -norma(1)
  norma(2)     = -norma(2)
  norma(ndime) = -norma(ndime)
 
  chkdi = dista
  
  coorj(1)     = xcoor(1)
  coorj(2)     = xcoor(2)
  coorj(ndime) = xcoor(ndime)


  do ipoib = 1,imbou(jimbo) % npoib       

     xcoor(1)     = imbou(jimbo) % cooib(1,ipoib) 
     xcoor(2)     = imbou(jimbo) % cooib(2,ipoib) 
     xcoor(ndime) = imbou(jimbo) % cooib(ndime,ipoib)
     
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
          imbou(jimbo) % npoib,mnoib,imbou(jimbo) % nboib,chkdi,&
          imbou(jimbo) % ltyib,imbou(jimbo) % lnoib,imbou(jimbo) % cooib,&
          tedis,tenor,proje,kboun,imbou(jimbo) % sabox,imbou(jimbo) % fabox,&
          imbou(jimbo) % blink,imbou(jimbo) % struc,imbou(jimbo) % ldist,&
          imbou(jimbo) % lnele)

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
  do iboib = 1,imbou(jimbo) % nboib
        do inoib = 1,mnoib
           ipoib = imbou(jimbo)%lnoib(inoib,iboib)      
                      
           bocod(1,inoib)     = imbou(jimbo) % cooib(1,ipoib) 
           bocod(2,inoib)     = imbou(jimbo) % cooib(2,ipoib) 
           bocod(ndime,inoib) = imbou(jimbo) % cooib(ndime,ipoib) 
        end do

        do inoib = 1,mnoib
           jnoib = inoib + 1
           if (jnoib > mnoib) jnoib = 1        
           
           call ibm_dsepar(iimbo,bocod(1,inoib),bocod(1,jnoib),chkdi,tedis,proji,projj,cooi1,cooi2)
           !
           ! A new shorest distance have been found (a tolerancy is used to prevent colinear segments)
           !
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
                 vect2(idime) = cooi2(idime) - cooi1(idime)
              end do
              call vecpro(vect1,vect2,norma,3)
              call vecuni(ndime,norma,dummi)
              
              tenor(1)     = 0.0_rp
              tenor(2)     = 0.0_rp
              tenor(ndime) = 0.0_rp
              do ilist =    1,imbou(jimbo)%lnele(  imbou(jimbo)%lnoib(inoib,iboib)  )%nelem
                 do jlist = 1,imbou(jimbo)%lnele(  imbou(jimbo)%lnoib(jnoib,iboib)  )%nelem
                    jboib =   imbou(jimbo)%lnele(  imbou(jimbo)%lnoib(inoib,iboib)  )%eleme(ilist)
                    kboib =   imbou(jimbo)%lnele(  imbou(jimbo)%lnoib(jnoib,iboib)  )%eleme(jlist)
                    if (jboib == kboib) then
                       call extbou(1_ip,nnode(imbou(jimbo)%ltyib(jboib)),imbou(jimbo) %lnoib(1,jboib),imbou(jimbo) %cooib,vect1)
                       do idime = 1,ndime  
                          tenor(idime) = tenor(idime) + vect1(idime)
                       end do
                       call vecuni(ndime,tenor,dummi)
                    end if
                 end  do
              end  do
              facto = 0.0_rp
              do idime = 1,ndime
                 facto = facto + norma(idime)*tenor(idime)
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
     print *,dista,iimbo,jimbo
     print *,"jopeta total en las partículas!!!!"
  !   pause
  end if

   
end subroutine ibm_dpapar



