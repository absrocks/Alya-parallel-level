require 'BOAST'
include BOAST
require_relative './KSplitOssBoast.rb'
require 'narray_ffi'
require '../Common/CommonArgs.rb'

class Example
  # To retrieve common arguments settings
  include CommonArgs    
  
  def self.run(vector_size)

    # Specify the kernel options
    #nests = [1,2,3,4,5,6,7,8,9,10,11,12]
		nests = (1..12).to_a
    #vector_size=4 #2 (C version)
    k_boast_params = {:vector_length => vector_size, :nests => nests, :unroll => false, :inline => :inline}
    
    set_lang(FORTRAN)
    #set_lang(C)

    # Create the kernel, ndime is needed to create the BOAST version
    k = KSplitBoast::new(k_boast_params)
    # Generate the will create a CKernel instance available using the variable kernel.
    # And can be used directly used as usuall 
    k.generate
    puts k.kernel   
    #k.kernel.build(:FCFLAGS => "-fimplicit-none -O2", :LDFLAGS => "-lgfortran")
		k.kernel.build(:FCFLAGS => "-fbounds-check -fimplicit-none", :LDFLAGS => "-lgfortran")
    #k.kernel.build(:FCFLAGS => "-cpp")
    
    # The CommonArgs module needs to be initialized.
    # We need to specify the number of kernel instances 
    # this way the module can create the correct number 
    # of arguments.
    ndime = 3
    seed = 10
    pgaus = 8
    pnode = 8
    CommonArgs.init(pgaus,pnode,ndime,vector_size,1,seed)
    
    # Arguments can still be set manually especially for those which the path 
    # taken depends of their value
    @@kfl_lumped = 1
    @@kfl_limit_nsi = 1
    @@kfl_stabi_nsi = 1

    # Then just run the kernel using the reference to the arguments
    #puts k.kernel.run(@@vsize, @@kfl_lumped, @@ndime,
    puts k.kernel.run(@@vsize, @@kfl_lumped, @@ndime,
       @@mnode, @@ntens, @@kfl_stabi_nsi,
       @@fvins_nsi, @@kfl_regim_nsi,
       @@kfl_press_nsi,@@kfl_linea_nsi, @@pabdf_nsi, @@nbdfp_nsi,
       @@kfl_sgsti_nsi, @@kfl_nota1_nsi, @@kfl_limit_nsi,
       @@penal_nsi, @@kfl_convection_type_nsi, @@nsi_galerkin,
       @@nsi_algebraic_split_oss, @@nsi_fractional_step_int,
       @@nsi_convection_conservative, @@nsi_convection_skew,@@nsi_convection_emac,
       @@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
       @@gpsha,@@gpcar,@@gpadv,@@gpvep[0],@@gpgrp[0],@@gprhs[0],@@gprhc[0],@@gpvel,
       @@gpgve,@@gpsgs,@@elvel,@@elpre,@@elbub,@@elauu[0],@@elaup[0],@@elapp[0],
       @@elapu[0],@@elrbu[0],@@elrbp[0],@@dtinv_loc,@@dtsgs,@@pbubl,
       @@gpsha_bub,@@gpcar_bub,@@elauq[0],@@elapq[0],@@elaqu[0],@@elaqp[0],@@elaqq[0],
       @@elrbq[0],@@densi)


  filename="nsi_element_assembly_split_oss_boast.f90"
	print filename
  File::open(filename, "w") { |f|
      f.puts k.kernel
  }

  puts "Done"
	return k

  end
end



options = {}

opt_parser = OptionParser.new { |opts|
  opts.banner = "Usage: run.rb --[options]=[value]"

  opts.on("-vVAL","--vector_size=VAL", "Specify vector_size value") { |n|
    options[:vector_size] = n.to_i
  }

  opts.on("-h", "--help", "Prints this help") {
    puts opts
    exit
  }
}.parse!

k = Example.run(options[:vector_size])
