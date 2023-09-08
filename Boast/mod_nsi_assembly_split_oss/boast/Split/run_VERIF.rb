require_relative './KSplitOssRef.rb'
require_relative './KSplitOssBoast.rb'
require 'narray_ffi'
require_relative '../Common/CommonArgs.rb'
require 'yaml'
require 'pp'
require 'csv'
require 'optparse'
require_relative '../Tools/LogInfo.rb'

class Experiment
  include CommonArgs

  def self.run(output_info)

    nests = (1..12).to_a
    #nests = (1..3).to_a
    
    # Recording platform informations
     LogInfo.init(output_info)
     LogInfo.get_info

    seed = 10
    epsilon = 10e-15
    vector_size=4 #2
		ndime = 3 
    charac = {:pgaus => 8, :pnode => 8}

    kernels = {}
    stats = {}


  #nests.each{|n|

      # This is the way we differenciate between ref and boast versions
      k_ref_params = {:kernel => :ref,:vector_length => vector_size, :preprocessor => false, :nests => nests, :unroll => false, :inline => :included, :CFLAGS => "-O2", :FFLAGS => "-O3"}
      k_boast_params = {:kernel => :boast, :vector_length => vector_size,:preprocessor => false, :nests => nests, :unroll => false, :inline=> :included, :CFLAGS => "-O2", :FFLAGS => "-O3", :FCFLAGS => "-cpp"}

      set_lang(FORTRAN)
      #set_lang(C)


			#.. REF KERNEL ..#        
      kernels[k_ref_params] = KSplitOssRef::new(k_ref_params)
      kernels[k_ref_params].generate
      kernels[k_ref_params].kernel.build(:FCFLAGS => "-fbounds-check -fimplicit-none", :LDFLAGS => "-lgfortran")
      #kernels[k_ref_params].kernel.build(:FCFLAGS => "-cpp" )
      LogInfo.register_kernel_info(k_ref_params, kernels[k_ref_params].kernel.to_s)

			#.. BOAST KERNEL ..#
      kernels[k_boast_params] = KSplitBoast::new(k_boast_params)
      kernels[k_boast_params].generate
      kernels[k_boast_params].kernel.build(:FCFLAGS => "-fbounds-check -fimplicit-none", :LDFLAGS => "-lgfortran")
      LogInfo.register_kernel_info(k_boast_params, kernels[k_boast_params].kernel.to_s)
     
 
      stats[k_ref_params] = {:time => []}
      stats[k_boast_params] = {:time => []}
      stats[k_ref_params][:characteristics] = charac
      stats[k_boast_params][:characteristics] = charac

      # Setting the values of the arguments
      CommonArgs.init(charac[:pgaus],charac[:pnode],ndime,vector_size,3,seed)


      @@kfl_lumped = 2 # 1
      @@kfl_limit_nsi = 1 # 2
      @@kfl_stabi_nsi = 1 # -1

		

      # Warm-up loop. Need to check if 100 is enough.
      2.times{|i|

				kernels[k_ref_params].kernel.run(@@vsize, @@kfl_lumped, @@ndime,
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


        kernels[k_boast_params].kernel.run(@@vsize, @@kfl_lumped, @@ndime,
       @@mnode, @@ntens, @@kfl_stabi_nsi,
       @@fvins_nsi, @@kfl_regim_nsi,
       @@kfl_press_nsi,@@kfl_linea_nsi, @@pabdf_nsi, @@nbdfp_nsi,
       @@kfl_sgsti_nsi, @@kfl_nota1_nsi, @@kfl_limit_nsi,
       @@penal_nsi, @@kfl_convection_type_nsi, @@nsi_galerkin,
       @@nsi_algebraic_split_oss, @@nsi_fractional_step_int,
       @@nsi_convection_conservative, @@nsi_convection_skew,@@nsi_convection_emac,
       @@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
       @@gpsha,@@gpcar,@@gpadv,@@gpvep[1],@@gpgrp[1],@@gprhs[1],@@gprhc[1],@@gpvel,
       @@gpgve,@@gpsgs,@@elvel,@@elpre,@@elbub,@@elauu[1],@@elaup[1],@@elapp[1],
       @@elapu[1],@@elrbu[1],@@elrbp[1],@@dtinv_loc,@@dtsgs,@@pbubl,
       @@gpsha_bub,@@gpcar_bub,@@elauq[1],@@elapq[1],@@elaqu[1],@@elaqp[1],@@elaqq[1],
       @@elrbq[1],@@densi)

     }
      

      2.times{|i|
        CommonArgs.init(charac[:pgaus],charac[:pnode],ndime,vector_size,3,seed)
        @@kfl_lumped = 2 # 1
        @@kfl_limit_nsi = 1 # 2
        @@kfl_stabi_nsi = 1 # -1
        
        stats[k_ref_params][:time][i] = kernels[k_ref_params].kernel.run(@@vsize, @@kfl_lumped, @@ndime,
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
       @@elrbq[0],@@densi)[:duration]


        stats[k_boast_params][:time][i] = kernels[k_boast_params].kernel.run(@@vsize, @@kfl_lumped, @@ndime,
       @@mnode, @@ntens, @@kfl_stabi_nsi,
       @@fvins_nsi, @@kfl_regim_nsi,
       @@kfl_press_nsi,@@kfl_linea_nsi, @@pabdf_nsi, @@nbdfp_nsi,
       @@kfl_sgsti_nsi, @@kfl_nota1_nsi, @@kfl_limit_nsi,
       @@penal_nsi, @@kfl_convection_type_nsi, @@nsi_galerkin,
       @@nsi_algebraic_split_oss, @@nsi_fractional_step_int,
       @@nsi_convection_conservative, @@nsi_convection_skew,@@nsi_convection_emac,
       @@pnode,@@pgaus,@@gpden,@@gpvis,@@gppor,@@gpsp1,@@gpsp2,@@gpvol,
       @@gpsha,@@gpcar,@@gpadv,@@gpvep[1],@@gpgrp[1],@@gprhs[1],@@gprhc[1],@@gpvel,
       @@gpgve,@@gpsgs,@@elvel,@@elpre,@@elbub,@@elauu[1],@@elaup[1],@@elapp[1],
       @@elapu[1],@@elrbu[1],@@elrbp[1],@@dtinv_loc,@@dtsgs,@@pbubl,
       @@gpsha_bub,@@gpcar_bub,@@elauq[1],@@elapq[1],@@elaqu[1],@@elaqp[1],@@elaqq[1],
       @@elrbq[1],@@densi)[:duration]


        diff_agrau = (@@agrau[0] - @@agrau[1]).abs
        diff_wgrgr = (@@wgrgr[0] - @@wgrgr[1]).abs
        diff_elauu = (@@elauu[0] - @@elauu[1]).abs
        diff_elrbu = (@@elrbu[0] - @@elrbu[1]).abs
        diff_elapu = (@@elapu[0] - @@elapu[1]).abs
        diff_elaqu = (@@elaqu[0] - @@elaqu[1]).abs
        diff_elaqp = (@@elaqp[0] - @@elaqp[1]).abs
        diff_elaqq = (@@elaqq[0] - @@elaqq[1]).abs
        diff_elapq = (@@elapq[0] - @@elapq[1]).abs
        diff_elauq = (@@elauq[0] - @@elauq[1]).abs
        diff_elaup = (@@elaup[0] - @@elaup[1]).abs
        diff_elapp = (@@elapp[0] - @@elapp[1]).abs
        diff_elrbp = (@@elrbp[0] - @@elrbp[1]).abs
        diff_elrbq = (@@elrbq[0] - @@elrbq[1]).abs

        diff_gpgrp = (@@gpgrp[1] - @@gpgrp[0]).abs
        diff_gprhs = (@@gprhs[1] - @@gprhs[0]).abs
        diff_gpvep = (@@gpvep[1] - @@gpvep[0]).abs
        diff_gprhc = (@@gprhc[1] - @@gprhc[0]).abs


        raise "Error: residue too big for agrau" if (diff_agrau > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for wgrgr" if (diff_wgrgr > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elauu" if (diff_elauu > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elrbu" if (diff_elrbu > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elapu" if (diff_elapu > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elaqu" if (diff_elaqu > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elaqp" if (diff_elaqp > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elaqq" if (diff_elaqq > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elapq" if (diff_elapq > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elauq" if (diff_elauq > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elaup" if (diff_elaup > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elapp" if (diff_elapp > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elrbp" if (diff_elrbp > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for elrbq" if (diff_elrbq > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for gpgrp" if (diff_gpgrp > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for gprhs" if (diff_gprhs > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for gprhc" if (diff_gprhc > epsilon).to_a.flatten.include? 1
        raise "Error: residue too big for gpvep" if (diff_gpvep > epsilon).to_a.flatten.include? 1
   # }
  }

    # Dumping info into a file
    LogInfo.dump_info
    puts "Done"
    return stats,kernels
  end
end

options = {}

opt_parser = OptionParser.new { |opts|
  opts.banner = "Usage: run.rb --[options]=[value]"

  opts.on("-dVAL", "--data=VAL", "Specify the path where to store the data in yaml format") { |n|
    options[:data_path] = n
  }

  opts.on("-kVAL", "--kernel=VAL", "Specify the path to store the kernel sources") { |n|
    options[:kernel_path] = n
  }

  opts.on("-iVAL", "--info=VAL", "Specify the path where to store the information such as the generated source, compiler info, etc... in yaml format") { |n|
    options[:info_path] = n
  }


  opts.on("-h", "--help", "Prints this help") {
    puts opts
    exit
  }
}.parse!

stats, kernels = Experiment.run(options[:kernel_path])

# In the results:
# Ref
# Boast

if options[:data_path] then
  File::open( options[:data_path], "w") { |f|
    f.print YAML::dump(stats)
  }
end

if options[:kernel_path] then
  File::open( options[:kernel_path], "w") { |f|
    kernels.each { |key, value|
      f.puts value.kernel if key[:kernel] == :boast
      #f.puts value.kernel if key[:kernel] == :ref
    }
  }
end
