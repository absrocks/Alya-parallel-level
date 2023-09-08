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

  def self.run(output_info, lang)

    nests = (1..12).to_a

    # Recording platform informations
     LogInfo.init(output_info)
     LogInfo.get_info

    seed = 10
    epsilon = 10e-15
    charac = {:pgaus => 8, :pnode => 8}
    repeat = 2

    kernels = {}
    stats = {}
    opt_space = OptimizationSpace::new( inline: [:inlined, :call, :included],
                                                unroll: [true, false],
                                                OFLAGS: ["-O2", "-O3"],
																								:vector_size => [2,4]
                                      )
    if lang == "C" then
      lang=C
    else
      lang=FORTRAN
    end

    set_lang(lang)


    optimizer = BruteForceOptimizer::new(opt_space, :randomize => true)


    puts optimizer.optimize { |opt|

      p opt

      k_ref_params = {:kernel => :ref, :vector_length => opt[:vector_size], :preprocessor => false, :nests => nests, :unroll => false, :inline => :included, :FCFLAGS => "#{opt[:OFLAGS]}", :CFLAGS => "#{opt[:OFLAGS]} -cpp", :LDFLAGS => "-lgfortran"}
      k_boast_params = {:kernel => :boast, :vector_length => opt[:vector_size], :preprocessor => false, :nests => nests, :unroll => opt[:unroll], :inline => opt[:inline], :FCFLAGS => "#{opt[:OFLAGS]}",:CFLAGS => "#{opt[:OFLAGS]} -cpp", :LDFLAGS => "-lgfortran"} 

      #.. REF KERNEL ..#      
      kernels[k_ref_params] = KSplitOssRef::new(k_ref_params)
      kernels[k_ref_params].generate
      kernels[k_ref_params].kernel.build(:LDFLAGS => k_ref_params[:LDFLAGS], :FCFLAGS => k_ref_params[:FCFLAGS])
      stats[k_ref_params]={:time => []}
      stats[k_ref_params][:characteristics] = charac
      LogInfo.register_kernel_info(k_ref_params, kernels[k_ref_params].kernel.to_s)


      #.. BOAST KERNEL ..#
      kernels[k_boast_params] = KSplitBoast::new(k_boast_params)
      kernels[k_boast_params].generate
      kernels[k_boast_params].kernel.build(:LDFLAGS => k_boast_params[:LDFLAGS], :FCFLAGS => k_boast_params[:FCFLAGS])
      stats[k_boast_params]={:time => []}
      stats[k_boast_params][:characteristics] = charac
      LogInfo.register_kernel_info(k_boast_params, kernels[k_boast_params].kernel.to_s)

			# Setting the values of the arguments
			ndime=3
      CommonArgs.init(charac[:pgaus],charac[:pnode],ndime,opt[:vector_size],3,seed)
      @@kfl_lumped = 2 # 1
      @@kfl_limit_nsi = 1 # 2
      @@kfl_stabi_nsi = 1 # -1


      # Warm-up loop. Need to check if 100 is enough.
      repeat.times{|i|
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
     }
     min_ref = stats[k_ref_params][:time].min
     p "Ref:  #{min_ref}"
		 min_ref


     repeat.times{|i|
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
     }
     min_boast = stats[k_boast_params][:time].min
     p "Boast:  #{min_boast}"
	   min_boast

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

  opts.on("-iVAL", "--info=VAL", "Specify the path where to store the information such as the generated source, compiler info, etc... in yaml format") { |n|
    options[:info_path] = n
  }

  opts.on("-lVAL","--lang=VAL", "Specify the language") { |n|
    options[:lang] = n
  }

  opts.on("-h", "--help", "Prints this help") {
    puts opts
    exit
  }
}.parse!


stats, kernels = Experiment.run(options[:info_path], options[:lang])
#stats, info = Experiment.run(options[:info_path],options[:dimension])

# In the results:
# Ref
# Boast

if options[:data_path] then
  File::open( options[:data_path], "w") { |f|
    f.print YAML::dump(stats)
  }
end
