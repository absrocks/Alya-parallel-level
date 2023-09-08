subroutine hlm_countobs_sites(delta,shot)

  use def_parame
  use def_master
  use def_domain
  use def_helmoz
  implicit none
  
  real(rp), intent(in) :: delta, shot

  integer(ip)         :: ielem,igaus,idime,inode,ii              !Indices and dimensions
  integer(ip)         :: iboun,inodb              !Indices and dimensions
  integer(ip)         :: pelty,pmate,pnode,ipoin,poin,bpoin
  integer(ip)         :: pgaus,plapl,porde,ptopo
  integer(ip)         :: pblty,pmatb,pnodb,pgaub

  real(rp)            :: accum
  integer(ip)         :: iobs
  real(rp),dimension(3,nsite_hlm) :: vobs
  complex(rp),dimension(nsite_hlm) :: dobs
  real(rp),dimension(nsite_hlm) :: sumobs 
  integer(ip)              :: iclos, iclosnode
  real(rp) :: normcoordsite,rclosdistance, rnodedistance 

  if(ISLAVE)then  

     do iobs=1,nsite_hlm 
        countobs(iobs)=0_ip
        sumobs(iobs)=0.0_rp
     end do

     elements: do ielem = 1,nelem
        pelty = ltype(ielem)       !Element type	
        pnode = nnode(pelty)       !Number of element nodes 
        pgaus = ngaus(pelty)       !Number of element Gauss points
        plapl = 0_ip               !Existence of Laplasian in formulation	
        porde = lorde(pelty)       !Element order
        ptopo = ltopo(pelty)       !Element topology
        pmate = 1_ip
        if ( nmate > 1_ip ) then
        pmate = lmate(ielem)
        end if
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do iobs=1,nsite_hlm
              iclosnode =clsite1_hlm(iobs) 
              if( lninv_loc(ipoin) ==  iclosnode ) then       
                 countobs(iobs) = countobs(iobs) + 1_ip
              end if
           end do
        end do
     end do elements

  else
     do iobs=1,nsite_hlm
        countobs(iobs)=0_ip
     end do
  end if

  call parari('SUM',0_ip,nsite_hlm,countobs)

end subroutine hlm_countobs_sites

