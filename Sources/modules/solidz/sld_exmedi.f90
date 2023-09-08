subroutine sld_exmedi(itask,pnode,pgaus,gpsha,gptlo,ielem_input)
  !------------------------------------------------------------------------
  !****f* Solidz/sld_exmedi
  ! NAME 
  !    sld_exmedi
  ! DESCRIPTION
  !    Coupling Exmedi-Solidz subroutine. 
  !    When the wave passes trough an element, the local time of reaction for each GP starts
  !    and is returned in GPTLO 
  !    
  !    FISOC     ... Isochrones vector (global-defined by EXMEDI)
  !    KACTI_SLD ... Activation vector (vector of "npoin" entries)
  !    TAULO     ... Local time of passage of the reaction at nodes.........tau 
  !    GPTLO     ... TAULO at gauss points..................................tau_GP
  ! 

  ! USES
  ! -
  ! USED BY
  ! sld_elmope
  !------------------------------------------------------------------------

  use def_master   ! general global variables
  use def_domain   ! geometry information
  use def_solidz   ! general solidz module information
  
  implicit none
  integer(ip)                :: jelem,igaus,kfiso,ivoig,ipoin,kpoin,inode,pelty,pnaux,pmate,jpoin,jnode
  integer(ip), intent(in)    :: pgaus,ielem_input,itask,pnode
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(out)   :: gptlo(mgaus)



  if (itask==0) then  !computing the activation vector
     
     ! reset activation vector
     do ipoin = 1,npoin
        kacti_sld(ipoin) = 0
     end do
 
     ! compute activation vector checking fisoc
     do jelem= 1,nelem
        pelty=ltype(jelem) 
        pnaux=nnode(pelty) 
        pmate = 1
        if( nmate_sld > 1 ) then
           pmate = lmate_sld(jelem)
        end if
        do inode=1,pnaux
           ipoin= lnods(inode,jelem) 
           if ((fisoc(ipoin) > 0.0_rp ) ) then 
              !
              ! fisoc is the cutim at which the wave passes a point
              ! fisoc is initialized to -1.0_rp
              !                               
              kacti_sld(ipoin) = kacti_sld(ipoin) + 1     !          
           end if                    
        end do        
     end do      
     !write(6,*) kacti_sld
     ! interchange kacti_sld (when running in parallel)
     call sld_parall(4)
   
     ! compute taulo (when running in parallel, this is done locally 
     ! because kacti_sld has already been interchanged)
     if (kfl_cellmod(pmate) == CELL_FITZHUGH_EXMEDI) then
        do ipoin = 1,npoin
           
           if (kacti_sld(ipoin) > 0) then 
              taulo(ipoin) = taulo(ipoin) + dtime                              
           end if
        end do
     end if
     
  else if (itask==1) then
 
     do igaus=1,pgaus
        gptlo(igaus) = 0.0_rp
        do inode=1,pnode
           ipoin= lnods(inode,ielem_input)
           gptlo(igaus)= gptlo(igaus) + taulo(ipoin)*gpsha(inode,igaus)         
        end do        
     end do
  end if

 


end subroutine sld_exmedi
