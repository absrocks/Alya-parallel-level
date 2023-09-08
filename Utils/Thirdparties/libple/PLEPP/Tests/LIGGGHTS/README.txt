#==========================================================================||==#
#===================================================| 2016Feb12@MARENOSTRUM |==#
+evolve 
  |_findCells 
    |_ locateM().findCell  


+ CFDEM/applications/solvers/cfdemSolverPiso/cfdemSolverPiso.C
  |_ CFDEM/src/lagrangian/cfdemParticle/cfdemCloud/cfdemCloud/cfdemCloud.C.particleCloud( mesh );
    |_ getDEMdata 
      |_ CFDEM/src/lagrangian/cfdemParticle/subModels/dataExchangeModel/twoWayMPI/twoWayMPI.C.dataExchangeM().getData  
        |_ data_liggghts_to_of  <- lmp 
        |_ data_of_to_liggghts  <- lmp  
      |_ CFDEM/src/lagrangian/cfdemParticle/subModels/dataExchangeModel/twoWayMPI/twoWayMPI.C.twoWayMPI
        |_ lmp = LAMMPS(...,comm_liggghts)
        |_ DEMts_ = lmp->update->dt;



+ data_liggghts_to_of | data_of_to_liggghts 
  |_ LIGGGHTS-PUBLIC/src/library_cfd_coupling.cpp.FixCfdCoupling.get_dc ->push | ->pull 
     |_ fix_cfd_coupling.h.get_dc -> CfdDatacoupling  
       |_ LIGGGHTS-PUBLIC/src/cfd_datacoupling_mpi.cpp.CfdDatacouplingMPI.push 
       |_ LIGGGHTS-PUBLIC/src/cfd_datacoupling_mpi.cpp.CfdDatacouplingMPI.pull
         |_ LIGGGHTS-PUBLIC/src/cfd_datacoupling_mpi.h.push_mpi  
           |_ from = find_push_property 
           |_ MPI_Allreduce  
         |_ LIGGGHTS-PUBLIC/src/cfd_datacoupling_mpi.h.pull_mpi
           |_ to = find_pull_property
           |_ MPI_Allreduce


#==========================================================================||==#
#===================================================| 2016Feb12@MARENOSTRUM |==#
+ main
  |_ LAMMPS *lammps = new LAMMPS(argc,argv,MPI_COMM_WORLD);
     LAMMPS::LAMMPS(...)
       memory   = new Memory(this);
       error    = new Error(this);
       universe = new Universe(this,communicator);

       MPI_Comm_split(universe->uworld,universe->iworld,0,&world);
       MPI_Comm_rank(world,&me);
       create()
         comm   = new Comm(this);
         domain = new Domain(this);
         force  = new Force(this);    // must be after group, to create temperature
         update = new Update(this);  // must be after output, force, neighbor
       init()
         update->init();
         force->init();         // pair must come after update due to minimizer
         domain->init();
         atom->init();          // atom must come after force and domain
         modify->init();        // modify must come after update, force, atom, domain
         neighbor->init();      // neighbor must come after force, modify
         comm->init();          // comm must come after force, modify, neighbor, atom
         output->init();        // output must come after domain, force, modify

+ Fix -> FixCfdCoupling
  |_ fix_cfd_coupling.h.FixCfdCoupling(lmp, ...)

FixCfdCoupling                  : public Fix
FixCfdCouplingConvection        : public Fix
FixCfdCouplingConvectionSpecies : public Fix
FixCfdCouplingForce             : public Fix


#==========================================================================||==#
#===================================================| 2016Feb15@MARENOSTRUM |==#

+ fix_cfd_coupling.cpp
 |_FixCfdCoupling::FixCfdCoupling
    |_ define CfdDataCouplingStyle(key,Class) else if (strcmp(arg[iarg_],#key) == 0) dc_ = new Class(lmp,iarg_+1,narg,arg,this);

+ FixCfdCoupling::end_of_step() "CFD Coupling established at step"  
+ CfdDatacouplingFile::readVectorData "Fix couple/cfd/file: Data corruption"


+ CfdDatacouplingFile::exchange
|_ CfdDatacouplingFile::push(pushnames_[i], pushtypes_[i], dummy, "");
|_ CfdDatacouplingFile::pull(pullnames_[i], pulltypes_[i], dummy, "");


+ CfdDatacouplingFile::pull
|_ to = CfdDatacoupling.find_pull_property
|_ readXXXData  ->  'dragforce', 'volumeweight'  

+ FixCfdCouplingForce::init()
|_ fix_coupling_->add_push_property #     //  values to be transfered to OF
|_ fix_coupling_->add_pull_property #     // values to come from OF
  
+ Universe::Universe

#==========================================================================||==#
#===================================================| 2016Feb17@MARENOSTRUM |==#

|_ lammps.cpp -> LAMMPS.init()
  |_ neighbor.cpp -> Neighbor.init() 
    |_ neighbor.cpp -> Neighbor.choose_build 
      |_ neigh_gran.cpp -> Neighbor.granular_bin_no_newton

|_ fix_cfd_coupling.cpp -> end_of_step 
  |_ CfdDatacoupling.exchange() 


|_ neighbor.h 
  |_ Neighbor::neigh_half_bin.cpp -> half_bin_newton_tri
    ibin = coord2bin(x[i]);
    ...
    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    list->inum = inum;

#==========================================================================||==#
#===================================================| 2016Feb19@MARENOSTRUM |==#
TODO:
  /home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016Feb15/Thirdparties/libple/PLEPP/Include/build_octree.c -> coord2bin
  /home/bsc21/bsc21704/z2016/REPOSITORY/LIGGGHTS/LIGGGHTS_2016FEB05/src/cfd_datacoupling_commdom.h -> set_bin_structure 

#==========================================================================||==#
#===================================================| 2016Feb24@MARENOSTRUM |==#
+ stopped!! :(  
|_ 1) LIGGGHTS_2016FEB05/src/cfd_datacoupling_commdom.h.set_bin_pts() because it is really slow!!
     a) it is necessary to use the coord2bin instead to use the "clasical" way to locate...   
     b) it is necessary the a meshless interpolation module 
     c) alya doesnt have a multiphase solver necessary to solve the fluid part...

+ Examples: 
  LIGGGHTS_2016FEB05/examples/LIGGGHTS/Tutorials_public/cohesion
  LIGGGHTS_2016FEB05/examples/LIGGGHTS/Tutorials_public/sph_2 

#==========================================================================||==#
#
#  LPP.x dump*
#

#==========================================================================||==#
#=====================================================| 2016ABR24@MINOTAURO |==#

FROM=/home/bsc21/bsc21704/z2016/REPOSITORY/LIGGGHTS/LIGGGHTS_2016FEB05/src
/home/bsc21/bsc21704/z2017/REPOSITORY/LIGGGHTS/LIGGGHTS_2016FEB05/examples/LIGGGHTS/Tutorials_public/cohesion



