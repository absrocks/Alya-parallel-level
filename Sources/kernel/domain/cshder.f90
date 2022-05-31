!------------------------------------------------------------------------
!> @addtogroup Domain 
!> @{
!> @file    cshder.f90
!> @author  Guillaume Houzeaux
!> @date    10/10/1972
!> @brief   Shape functions, derivatives and Hessians
!> @details Shape functions, derivatives and Hessians:
!>          \verbatim
!>          - For each element type, using user integration rule:
!>            WEIGP(ngaus)
!>            SHAPE(nnode,ngaus)
!>            DERIV(ndime,nnode,ngaus)
!>            HESLO(ntens,nnode,ngaus)
!>          - For each element type, the bubble:
!>            SHAPE_BUB(ngaus)
!>            DERIV_BUB(ndime,ngaus)
!>            HESLO_BUB(ntens,ngaus)
!>          - For each element type, using a closed integration rule:
!>            WEIGC(nnode)
!>            SHAPC(nnode,nnode)
!>            DERIC(ndime,nnode,nnode)
!>            HESLC(ntens,nnode,nnode)
!>          - Center of gravity:
!>            SHACG(nnode)
!>            DERCG(ndime,nnode)
!>            WEICG
!>          - Element Gauss points to nodes: 
!>            SHAGA(ngaus,nnode) 
!>          \endverbatim
!> @} 
!------------------------------------------------------------------------

subroutine cshder(itask)
  use def_master, only : kfl_paral
  use def_kintyp, only : ip,rp
  use def_domain, only : elmar
  use def_domain, only : nnode
  use def_domain, only : ngaus
  use def_domain, only : ldime
  use def_domain, only : lrule
  use def_domain, only : lquad
  use def_domain, only : ltopo
  use def_domain, only : lexis
  use def_domain, only : mnode
  use def_domain, only : nelty
  use def_domain, only : ntens
  use def_domain, only : ndime
  use def_domain, only : memor_dom
  use def_domain, only : nelem
  use def_domain, only : mgaus
  use def_domain, only : ltype
  use def_kermod, only : kfl_data_base_array
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use mod_elmgeo, only : elmgeo_shapf_deriv_heslo_bubble
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: igaus,ielty,inode,mdime,pelty,ielem
  integer(ip)             :: pnode,pgaus,pdime,prule,ptopo,pquad
  real(rp)                :: poscg(ndime) 
  real(rp)                :: posgw(ndime,mnode) 
  real(rp)                :: weigw(mnode)
  real(rp)                :: posgc(ndime,mnode)
  integer(ip)             :: ierro

  if( itask == 1 ) then

     !----------------------------------------------------------------------
     !
     ! Allocate memory for structure ELMAR
     !
     !----------------------------------------------------------------------

     allocate( elmar(nelty) )     
     do ielty = 1,nelty
        nullify( elmar(ielty) % shape     )
        nullify( elmar(ielty) % deriv     )
        nullify( elmar(ielty) % heslo     )
        nullify( elmar(ielty) % posgp     )    
        nullify( elmar(ielty) % weigp     )
        nullify( elmar(ielty) % shape_bub )
        nullify( elmar(ielty) % deriv_bub )
        nullify( elmar(ielty) % heslo_bub )
        nullify( elmar(ielty) % shacg     )
        nullify( elmar(ielty) % dercg     )
        nullify( elmar(ielty) % hescg     )   
        nullify( elmar(ielty) % shaga     )
        nullify( elmar(ielty) % shapc     )
        nullify( elmar(ielty) % deric     )
        nullify( elmar(ielty) % heslc     )
        nullify( elmar(ielty) % weigc     )    
        nullify( elmar(ielty) % shaib     )
        nullify( elmar(ielty) % derib     )
        nullify( elmar(ielty) % weiib     )
     end do

  else if( itask == 2 ) then 

     !----------------------------------------------------------------------
     !
     ! Deallocate memory for structure ELMAR
     !
     !----------------------------------------------------------------------

     do ielty = 1,nelty
        call memory_deallo(memor_dom,'SHAPE'    ,'cshder',elmar(ielty) % shape)
        call memory_deallo(memor_dom,'DERIV'    ,'cshder',elmar(ielty) % deriv)
        call memory_deallo(memor_dom,'HESLO'    ,'cshder',elmar(ielty) % heslo)
        call memory_deallo(memor_dom,'POSGP'    ,'cshder',elmar(ielty) % posgp)
        call memory_deallo(memor_dom,'WEIGP'    ,'cshder',elmar(ielty) % weigp)

        call memory_deallo(memor_dom,'SHAPE_BUB','cshder',elmar(ielty) % shape_bub)
        call memory_deallo(memor_dom,'DERIV_BUB','cshder',elmar(ielty) % deriv_bub)
        call memory_deallo(memor_dom,'HESLO_BUB','cshder',elmar(ielty) % heslo_bub)

        call memory_deallo(memor_dom,'SHACG'    ,'cshder',elmar(ielty) % shacg)
        call memory_deallo(memor_dom,'DERCG'    ,'cshder',elmar(ielty) % dercg)
        call memory_deallo(memor_dom,'HESCG'    ,'cshder',elmar(ielty) % hescg)

        call memory_deallo(memor_dom,'SHAGA'    ,'cshder',elmar(ielty) % shaga)

        call memory_deallo(memor_dom,'SHAPC'    ,'cshder',elmar(ielty) % shapc)
        call memory_deallo(memor_dom,'DERIC'    ,'cshder',elmar(ielty) % deric)
        call memory_deallo(memor_dom,'HESLC'    ,'cshder',elmar(ielty) % heslc)
        call memory_deallo(memor_dom,'WEIGC'    ,'cshder',elmar(ielty) % weigc)
     end do

     deallocate(elmar)

  else if( itask == 3 ) then

     !----------------------------------------------------------------------
     !
     ! Allocate: 3. Deallocate: 4
     !
     !----------------------------------------------------------------------

     ierro = 0
     poscg = 0.0_rp
     posgw = 0.0_rp
     weigw = 0.0_rp
     posgc = 0.0_rp
     !
     ! Element shape function and derivatives: SHAPE,DERIV,HESLO,WEIGP 
     !
     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then 
           
           pnode = nnode(ielty)
           pgaus = ngaus(ielty)
           pdime = ldime(ielty)
           mdime = max(1_ip,pdime)
           prule = lrule(ielty)
           elmar(ielty) % pgaus = pgaus
           
           call memory_alloca(memor_dom,'SHAPE'    ,'cshder',elmar(ielty) % shape,pnode,pgaus)
           call memory_alloca(memor_dom,'DERIV'    ,'cshder',elmar(ielty) % deriv,mdime,pnode,pgaus)
           call memory_alloca(memor_dom,'HESLO'    ,'cshder',elmar(ielty) % heslo,ntens,pnode,pgaus)
           call memory_alloca(memor_dom,'POSGP'    ,'cshder',elmar(ielty) % posgp,ndime,pgaus)
           call memory_alloca(memor_dom,'WEIGP'    ,'cshder',elmar(ielty) % weigp,pgaus)

           call memory_alloca(memor_dom,'SHAPE_BUB','cshder',elmar(ielty) % shape_bub,pgaus)
           call memory_alloca(memor_dom,'DERIV_BUB','cshder',elmar(ielty) % deriv_bub,mdime,pgaus)
           call memory_alloca(memor_dom,'HESLO_BUB','cshder',elmar(ielty) % heslo_bub,ntens,pgaus)

           call rulepw(pdime,pgaus,prule,elmar(ielty) % posgp,elmar(ielty) % weigp,ierro)
           call shafal(&
                elmar(ielty) % posgp,pdime,pnode,pgaus,ntens,elmar(ielty) % shape,elmar(ielty) % deriv,&
                elmar(ielty) % heslo,ierro)
      
           do igaus = 1,pgaus
              call elmgeo_shapf_deriv_heslo_bubble(&
                   mdime,pnode,elmar(ielty) % posgp(1:mdime,igaus),elmar(ielty) % shape_bub(igaus:igaus),&
                   elmar(ielty) % deriv_bub(1:mdime,igaus),elmar(ielty) % heslo_bub(1:ntens,igaus))
           end do

        end if
     end do
     !
     ! Element Center of gravity shape function and derivatives SHACG,DERCG,HESCG,WEICG 
     !
     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then

           pnode = nnode(ielty)
           pgaus = ngaus(ielty)
           pdime = ldime(ielty)
           mdime = max(1_ip,pdime)
           prule = lrule(ielty)

           call memory_alloca(memor_dom,'SHACG','cshder',elmar(ielty) % shacg,pnode)
           call memory_alloca(memor_dom,'DERCG','cshder',elmar(ielty) % dercg,mdime,pnode)
           call memory_alloca(memor_dom,'HESCG','cshder',elmar(ielty) % hescg,ntens,pnode)

           prule = ltopo(ielty) * 2 + 1 
           call rulepw(pdime,1_ip,prule,poscg,elmar(ielty)%weicg,ierro)
           call shafun(&
                poscg,pdime,pnode,ntens,&
                elmar(ielty)%shacg,elmar(ielty)%dercg,elmar(ielty)%hescg,ierro)

        end if
     end do
     !
     ! Computes the interpolation functions associated to the integration
     ! points in order to extrapolate values from these integration points
     ! to the nodes
     !
     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then

           pnode = nnode(ielty)
           pgaus = ngaus(ielty)
           pdime = ldime(ielty)
           pquad = lquad(ielty)
           mdime = max(1_ip,pdime)

           call memory_alloca(memor_dom,'SHAGA','cshder',elmar(ielty) % shaga,pgaus,pnode)

           if( lquad(ielty) == 1 ) then
              do inode = 1,pnode
                 do igaus = 1,pgaus
                    elmar(ielty) % shaga(igaus,inode) = 0.0_rp
                 end do
                 elmar(ielty) % shaga(inode,inode) = 1.0_rp
              end do
           else
              
              ptopo = ltopo(ielty) 
              prule = (ptopo+1)*2
              call rulepw(pdime,pnode,prule,posgw,weigw,ierro)
              call shafga(&
                   posgw,pdime,ptopo,pgaus,pnode,&
                   elmar(ielty)%shaga,ierro)
           end if

        end if
     end do
     !
     ! Element shape function and derivatives SHAPC,DERIC,HESLC,WEIGC 
     ! for a close rule
     !
     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then

           pnode = nnode(ielty)
           pdime = ldime(ielty)
           mdime = max(1_ip,pdime)

           call memory_alloca(memor_dom,'SHAPC','cshder',elmar(ielty) % shapc,pnode,pnode)
           call memory_alloca(memor_dom,'DERIC','cshder',elmar(ielty) % deric,mdime,pnode,pnode)
           call memory_alloca(memor_dom,'HESLC','cshder',elmar(ielty) % heslc,ntens,pnode,pnode)
           call memory_alloca(memor_dom,'WEIGC','cshder',elmar(ielty) % weigc,pnode)

           ptopo =ltopo(ielty)
           prule = (ptopo+1) * 2

           call rulepw(pdime,pnode,prule,posgc,elmar(ielty)%weigc,ierro)
           call shafal(&
                posgc,pdime,pnode,pnode,ntens,elmar(ielty)%shapc,&
                elmar(ielty)%deric,elmar(ielty)%heslc,ierro)

        end if
     end do

  end if

end subroutine cshder
