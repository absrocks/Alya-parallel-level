
require_relative '../Common/Parameters.rb'
class KSplitOss
  include Parameters
  attr_reader :kernel

  def declare_locs()
   Parameters.initialize(@opts[:vector_length])
   @locs = []
   @locs.push @inode	= $inode
   @locs.push @c1	= $c1
   @locs.push @c2	= $c2
   @locs.push @c3	= $c3
   @locs.push @c4	= $c4
   @locs.push @alpha	= $alpha
   @locs.push @beta	= $beta
   @locs.push @fact1_p	= $fact1_p
   @locs.push @jnode	= $jnode
   @locs.push @jdime	= $jdime
   @locs.push @idofv	= $idofv
   @locs.push @jdofv	= $jdofv
   @locs.push @ivect	= $ivect
   @locs.push @igaus	= $igaus
   @locs.push @idime	= $idime
   @locs.push @jdof	= $jdof
   @locs.push @itime	= $itime
   @locs.push @fact	= $fact
   @locs.push @idof	= $idof
   @locs.push @idof1	= $idof1
   @locs.push @jdof1	= $jdof1
   @locs.push @gpsp1_p	= $gpsp1_p
   @locs.push @gpsp1_v	= $gpsp1_v
   @locs.push @gpsp2_v	= $gpsp2_v
   @locs.push @dtinv_mod	= $dtinv_mod
   @locs.push @gpveo	= $gpveo
   #@locs.push @nsi_galerkin	= $nsi_galerkin
   #@locs.push @nsi_algebraic_split_oss	= $nsi_algebraic_split_oss
   #@locs.push @nsi_fractional_step	= $nsi_fractional_step
   @locs.push @wgrgr     = $wgrgr
   @locs.push @agrau     = $agrau
  end


  def declare_parameters()
    Parameters.initialize(@opts[:vector_length])
    @args = []

		
    @args.push @vsize    = $vsize
    @args.push @kfl_lumped    = $kfl_lumped
    @args.push @ndime       = $ndime

    #if ndime.nil? then
    #  @args.push @ndime       = $ndime
    #else
    #  @ndime                  = ndime
    #end

    @args.push @mnode     = $mnode
    @args.push @ntens     = $ntens
    @args.push @kfl_stabi_nsi = $kfl_stabi_nsi
    @args.push @fvins_nsi     = $fvins_nsi
    @args.push @kfl_regim_nsi = $kfl_regim_nsi
    @args.push @kfl_press_nsi = $kfl_press_nsi
    @args.push @kfl_linea_nsi = $kfl_linea_nsi
    @args.push @pabdf_nsi     = $pabdf_nsi
    @args.push @nbdfp_nsi     = $nbdfp_nsi
    @args.push @kfl_sgsti_nsi = $kfl_sgsti_nsi
    @args.push @kfl_nota1_nsi = $kfl_nota1_nsi
    @args.push @kfl_limit_nsi = $kfl_limit_nsi
    @args.push @penal_nsi     = $penal_nsi
    @args.push @kfl_convection_type_nsi     = $kfl_convection_type_nsi
    @args.push @nsi_galerkin     = $nsi_galerkin
    @args.push @nsi_algebraic_split_oss     = $nsi_algebraic_split_oss
    @args.push @nsi_fractional_step_int     = $nsi_fractional_step_int
    @args.push @nsi_convection_conservative     = $nsi_convection_conservative
    @args.push @nsi_convection_skew     = $nsi_convection_skew
    @args.push @nsi_convection_emac     = $nsi_convection_emac
    @args.push @pnode     = $pnode
    @args.push @pgaus     = $pgaus
    @args.push @gpden     = $gpden
    @args.push @gpvis     = $gpvis
    @args.push @gppor     = $gppor
    @args.push @gpsp1     = $gpsp1
    @args.push @gpsp2     = $gpsp2
    @args.push @gpvol     = $gpvol
    @args.push @gpsha     = $gpsha
    @args.push @gpcar     = $gpcar
    @args.push @gpadv     = $gpadv
    @args.push @gpvep     = $gpvep
    @args.push @gpgrp     = $gpgrp
    @args.push @gprhs     = $gprhs
    @args.push @gprhc     = $gprhc
    @args.push @gpvel     = $gpvel
    @args.push @gpgve     = $gpgve
    @args.push @gpsgs     = $gpsgs
    @args.push @elvel     = $elvel
    @args.push @elpre     = $elpre
    @args.push @elbub     = $elbub
    @args.push @elauu     = $elauu
    @args.push @elaup     = $elaup
    @args.push @elapp     = $elapp
    @args.push @elapu     = $elapu
    @args.push @elrbu     = $elrbu
    @args.push @elrbp     = $elrbp
    @args.push @dtinv_loc = $dtinv_loc
    @args.push @dtsgs     = $dtsgs
    @args.push @pbubl     = $pbubl
    @args.push @gpsha_bub = $gpsha_bub
    @args.push @gpcar_bub = $gpcar_bub
    @args.push @elauq     = $elauq
    @args.push @elapq     = $elapq
    @args.push @elaqu     = $elaqu
    @args.push @elaqp     = $elaqp
    @args.push @elaqq     = $elaqq
    @args.push @elrbq     = $elrbq
    @args.push @densi	  = $densi




    #@args.push @kfl_duatss    = $kfl_duatss
    #@args.push @fact_duatss   = $fact_duatss
    #@args.push @fcons_nsi     = $fcons_nsi
    #@args.push @bemol_nsi     = $bemol_nsi
    #@args.push @fvela_nsi     = $fvela_nsi
    #@args.push @kfl_rmom2_nsi = $kfl_rmom2_nsi
    #@args.push @kfl_p1ve2_nsi = $kfl_p1ve2_nsi
    #@args.push @kfl_confi_nsi = $kfl_confi_nsi
    #@args.push @kfl_penal_nsi = $kfl_penal_nsi
    #@args.push @kfl_bubbl_nsi = $kfl_bubbl_nsi


  end

  def initialize(options)
    @opts = {:vector_length => 1, :preprocessor => false, :nests => (1..12).to_a, :unroll => false, :inline => :included}
    @opts.update(options)
    declare_parameters()
    declare_locs()
  end

    def declare_procedure(name, functions = nil)
      return Procedure(name, @args, :functions => functions, :locals => @locs)
    end
end
