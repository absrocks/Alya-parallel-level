subroutine dod_outcpu
  !-----------------------------------------------------------------------
  !****f* dodeme/dod_outcpu
  ! NAME 
  !    dod_outcpu
  ! DESCRIPTION
  !    This routine outputs Dodeme CPU time table
  ! USES
  !
  ! USED BY
  !
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_dodeme
  use mod_outfor, only : outfor
  implicit none
  real(rp) :: onvto
  !
  ! Title
  !
  routp(1) = cpu_servi(1,servi)
  call outfor(29_ip,lun_outpu_dod,' ')
  if( cpu_servi(1,servi) <= 0.0_rp ) cpu_servi(1,servi) = 1.0_rp
  onvto = 100.0_rp / cpu_servi(1,servi)

  coutp(1) = 'COPY MESH TO SUBD. MESHES'
  routp(2) = cpu_dod_copy_mesh
  routp(3) = routp(2) * onvto
  call outfor(30_ip,lun_outpu_dod,' ')

  coutp(1) = 'SUBDOMAIN MESH GRAPHS'
  routp(2) = cpu_dod_graphs
  routp(3) = routp(2) * onvto
  call outfor(30_ip,lun_outpu_dod,' ')
 
  coutp(1) = 'SUBDOMAIN SKD-TREES'
  routp(2) = cpu_dod_kdtree
  routp(3) = routp(2) * onvto
  call outfor(30_ip,lun_outpu_dod,' ')

  coutp(1) = 'SUBDOMAIN ELEMENT SEARCH'
  routp(2) = cpu_dod_elsest
  routp(3) = routp(2) * onvto
  call outfor(30_ip,lun_outpu_dod,' ')
  !
  ! => HOLE CUTTING
  !
  coutp(1) = 'HOLE CUTTING'
  routp(2) = cpu_dod_holcut
  routp(3) = routp(2) * onvto
  call outfor(30_ip,lun_outpu_dod,' ')

  coutp(1) = 'MARK NODES'
  routp(2) = cpu_dod_holcut_marknodes
  routp(3) = routp(2) * onvto
  call outfor(31_ip,lun_outpu_dod,' ')

  coutp(1) = 'INVERSE HOLES'
  routp(2) = cpu_dod_holcut_inversehole
  routp(3) = routp(2) * onvto
  call outfor(31_ip,lun_outpu_dod,' ')

  coutp(1) = 'MARK HOLE ELEMENTS'
  routp(2) = cpu_dod_holcut_markelements
  routp(3) = routp(2) * onvto
  call outfor(31_ip,lun_outpu_dod,' ')

  coutp(1) = 'HOLES BOUNDARIES'
  routp(2) = cpu_dod_holcut_holeboundary
  routp(3) = routp(2) * onvto
  call outfor(31_ip,lun_outpu_dod,' ')

  coutp(1) = 'HOLES FRINGE NODES'
  routp(2) = cpu_dod_holcut_fringenodes
  routp(3) = routp(2) * onvto
  call outfor(31_ip,lun_outpu_dod,' ')

  coutp(1) = 'GRAPHS & SKD-TREES'
  routp(2) = cpu_dod_graphs_kdtree
  routp(3) = routp(2) * onvto
  call outfor(31_ip,lun_outpu_dod,' ')
  !
  ! => END HOLE CUTTING
  !
  coutp(1) = 'PATCH INTERFACES'
  routp(2) = cpu_dod_patch_boundaries
  routp(3) = routp(2) * onvto
  call outfor(30_ip,lun_outpu_dod,' ')

  coutp(1) = 'PRESCRIBED INTERFACES'
  routp(2) = cpu_dod_prescribed_boundaries
  routp(3) = routp(2) * onvto
  call outfor(30_ip,lun_outpu_dod,' ')

  coutp(1) = 'CANDIDATE NODES'
  routp(2) = cpu_dod_candidate_nodes
  routp(3) = routp(2) * onvto
  call outfor(30_ip,lun_outpu_dod,' ')

  coutp(1) = 'ENXTENSION'
  routp(2) = cpu_dod_extension
  routp(3) = routp(2) * onvto
  call outfor(30_ip,lun_outpu_dod,' ')

end subroutine dod_outcpu
