
require_relative '../Common/subroutine.rb'
require_relative './KSplitOss.rb'


class Nest1 < Subroutine
  def initialize(options,functions = nil)
    super("nest1",options,functions)

    # ins
    @mnode = Parameters.copy($mnode,:in)
    @pnode = Parameters.copy($pnode,:in)
    @pgaus = Parameters.copy($pgaus,:in)
    @gpden = Parameters.copy($gpden,:in)
    @gpcar = Parameters.copy($gpcar,:in)
    @gpadv = Parameters.copy($gpadv,:in)
    @ndime = Parameters.copy($ndime,:in)
    @vsize = Parameters.copy($vsize,:in)
    # inouts
    @wgrgr = Parameters.copy($wgrgr,:inout)
    @agrau = Parameters.copy($agrau,:inout)
		# ADD TO ARGS
    @args.push @vsize,@ndime,@mnode,@pnode,@pgaus,@gpden,@gpcar,@gpadv,@wgrgr,@agrau

    # locals
    @inode = Parameters.copy($inode)
    @jnode = Parameters.copy($jnode)
    @igaus = Parameters.copy($igaus)
    @idime = Parameters.copy($idime)
		# ADD TO LOCS
    @locs.push @inode, @jnode, @igaus, @idime
  end
 
  def generate


    block = lambda{
      
	
      pr For(@igaus,1,@pgaus){
        pr For(@inode,1,@pnode){

#						if get_lang == FORTRAN then
#      				pr @agrau[@inode,@igaus]===Set(0.0,@agrau[@inode,@igaus])
#      				pr @wgrgr[@inode,@inode,@igaus]===Set(0.0,@wgrgr[@inode,@inode,@igaus])
#            elsif get_lang == C
#              code =<<EOF
#              memset(&#{@agrau[@inode,@igaus]},0,sizeof(#{@agrau[@inode,@igaus]}));   
#              memset(&#{@wgrgr[@inode,@inode,@igaus]},0,sizeof(#{@wgrgr[@inode,@inode,@igaus]}));   
#EOF
#              get_output.print code
#            end

          pr For(@idime,1,@ndime){
            pr @agrau[@inode,@igaus] === @agrau[@inode,@igaus] + @gpadv[@idime,@igaus]*@gpcar[@idime,@inode,@igaus]
          }
          pr @agrau[@inode,@igaus] === @gpden[@igaus] * @agrau[@inode,@igaus]
          pr For(@jnode,1,@pnode){   

						if get_lang == FORTRAN then
            	pr @wgrgr[@inode,@jnode,@igaus] === 0.0
						elsif get_lang == C
              code =<<EOF
							memset(&#{@wgrgr[@inode,@jnode,@igaus]},0,sizeof(#{@wgrgr[@inode,@jnode,@igaus]})); 
EOF
							get_output.print code
            end

            pr For(@idime,1,@ndime){
              pr @wgrgr[@inode,@jnode,@igaus] === @wgrgr[@inode,@jnode,@igaus] + @gpcar[@idime,@inode,@igaus]*@gpcar[@idime,@jnode,@igaus]
            }
          }
        }
      }
    }
    construct(block)
  end
end

class Nest2 < Subroutine
  def initialize(options, functions = nil)
    super("nest2",options,functions)

    # ins
    @mnode = Parameters.copy($mnode,:in)
    @ndime = Parameters.copy($ndime,:in)
    @pnode = Parameters.copy($pnode,:in)
    @pgaus = Parameters.copy($pgaus,:in)
    @gpden = Parameters.copy($gpden,:in)
    @gpcar = Parameters.copy($gpcar,:in)
    @gpsp2 = Parameters.copy($gpsp2,:in)
    @gpvol = Parameters.copy($gpvol,:in)
    @gpvis = Parameters.copy($gpvis,:in)
    @dtinv_mod = Parameters.copy($dtinv_mod,:in)
    @gppor = Parameters.copy($gppor,:in)
    @gpsha = Parameters.copy($gpsha,:in)
    @gpsp1_v = Parameters.copy($gpsp1_v,:in)
    @wgrgr = Parameters.copy($wgrgr,:in)
    @agrau = Parameters.copy($agrau,:in)
    @vsize = Parameters.copy($vsize,:in)
    #inouts
    @elauu = Parameters.copy($elauu,:inout)
    @pabdf_nsi = Parameters.copy($pabdf_nsi, :inout)

    @args.push @vsize,@ndime,@mnode,@pnode,@pgaus,@gpden,@gpcar,@gpsp2,@gpvol,@gpvis,@dtinv_mod,
                @gppor,@gpsha,@gpsp1_v,@wgrgr,@agrau,@elauu,@pabdf_nsi
  end
  def generate

    # locals
    @fact = Parameters.copy($fact)
    @igaus = Parameters.copy($igaus)
    @inode = Parameters.copy($inode)
    @jnode = Parameters.copy($jnode)
    @idofv = Parameters.copy($idofv)
    @jdofv = Parameters.copy($jdofv)
    @idime = Parameters.copy($idime)
    @jdime = Parameters.copy($jdime)

    @locs.push @fact, @idofv, @jdofv, @igaus, @inode, @jnode, @idime, @jdime 

    for1 = For(@igaus,1,@pgaus){ 
      pr @fact[1] === @gpsp2[@igaus] * @gpvol[@igaus]
      pr @fact[7] === @gpvis[@igaus] * @gpvol[@igaus]
      pr @fact[8] === @gpsp1_v[@igaus] * @gpvol[@igaus]
      pr @fact[9] ===  @gpden[@igaus] * @pabdf_nsi[1] * @dtinv_mod + @gppor[@igaus]
    
      pr For(@inode,1,@pnode){
        pr For(@idime,1,@ndime){
          pr @idofv[1] === (@inode-1)*@ndime+@idime
          pr @fact[2] === @fact[1] * @gpcar[@idime,@inode,@igaus]      
          pr For(@jnode,1,@pnode){
            pr For(@jdime,1,@ndime){
              pr @jdofv[@jdime]              === (@jnode-1)*@ndime+@jdime
              pr @elauu[@idofv[1],@jdofv[@jdime]] === @elauu[@idofv[1],@jdofv[@jdime]] + @fact[2] * @gpcar[@jdime,@jnode,@igaus]                   
            }.unroll(@unroll)
    
            pr @jdofv[1]              === (@jnode-1)*@ndime+@idime
            pr @fact[5]              === @gpsha[@inode,@igaus] * @gpvol[@igaus]
            pr @fact[6]              === @fact[5] * ( @agrau[@jnode,@igaus] + @fact[9] * @gpsha[@jnode,@igaus] ) + @fact[7] * @wgrgr[@inode,@jnode,@igaus] + @fact[8] *   @agrau[@jnode,@igaus] * @agrau[@inode,@igaus]
            pr @elauu[@idofv[1],@jdofv[1]] === @elauu[@idofv[1],@jdofv[1]] + @fact[6]
    
          }
        }
      }
    }
    
    main_block = lambda{
      pr for1
    }
          

    construct(main_block)
  end
end




class Nest3 < Subroutine
  def initialize(options, functions = nil)
     super("nest3",options,functions)

     # ins
     @mnode = Parameters.copy($mnode,:in) 
     @pgaus = Parameters.copy($pgaus,:in)
     @pnode = Parameters.copy($pnode,:in)
     @gpden = Parameters.copy($gpden,:in)
     @ndime = Parameters.copy($ndime,:in)
     @gpvol = Parameters.copy($gpvol,:in)
     @gpgve = Parameters.copy($gpgve,:in)
     @gpvel = Parameters.copy($gpvel,:in)
     @gpadv = Parameters.copy($gpadv,:in)
     @gpcar = Parameters.copy($gpcar,:in)
     @gpsha = Parameters.copy($gpsha,:in)	
     @vsize = Parameters.copy($vsize,:in)
		
     @kfl_regim_nsi = Parameters.copy($kfl_regim_nsi,:in)
     @elvel = Parameters.copy($elvel,:in)
     @nbdfp_nsi = Parameters.copy($nbdfp_nsi,:in)
     @densi = Parameters.copy($densi,:in)
     @dtinv_loc = Parameters.copy($dtinv_loc,:in)
     @nsi_convection_skew = Parameters.copy($nsi_convection_skew,:in)
     @nsi_convection_emac = Parameters.copy($nsi_convection_emac,:in)
     @nsi_convection_conservative = Parameters.copy($nsi_convection_conservative,:in)
     @kfl_convection_type_nsi = Parameters.copy($kfl_convection_type_nsi,:in)
     # inouts
     @elauu = Parameters.copy($elauu,:inout)
     @elrbu = Parameters.copy($elrbu,:inout)
     @gpveo = Parameters.copy($gpveo,:inout)
     @pabdf_nsi = Parameters.copy($pabdf_nsi, :inout)

     @args.push @vsize,@mnode,@pgaus,@pnode,@gpden,@gpvol,@gpadv,@gpcar,@gpsha,@kfl_regim_nsi,@elvel,@nbdfp_nsi,@densi,@dtinv_loc,@elauu,@elrbu,@gpveo,@pabdf_nsi, @kfl_convection_type_nsi, @nsi_convection_conservative, @gpgve, @nsi_convection_emac,@ndime, @gpvel, @nsi_convection_skew
  end
  def generate
     # local
     @igaus = Parameters.copy($igaus)
     @inode = Parameters.copy($inode)
     @jnode = Parameters.copy($jnode)
     @idime = Parameters.copy($idime)
     @jdime = Parameters.copy($jdime)
     @fact = Parameters.copy($fact)
     @itime = Parameters.copy($itime)
     @idofv = Parameters.copy($idofv)
     @jdofv = Parameters.copy($jdofv)

     @locs.push @igaus, @inode, @jnode, @idime, @jdime, @itime, @fact, @idofv, @jdofv


     for_inode = For(@inode,1,@pnode)
     for_jnode = For(@jnode,1,@pnode)
     for_idime = For(@idime,1,@ndime)
     for_jdime = For(@jdime,1,@ndime)
     for_itime = For(@itime,2,@nbdfp_nsi)
		 for_igaus = For(@igaus,1,@pgaus)


		@cond1 = @kfl_convection_type_nsi == @nsi_convection_conservative
		@cond2 = @kfl_convection_type_nsi == @nsi_convection_skew
		@cond3 = @kfl_convection_type_nsi == @nsi_convection_emac


		@block1 = lambda {
					puts "CONSERVATIVE FORM NOT CODED"
					opn for_igaus
						pr @fact[1] === @gpden[@igaus] * @gpvol[@igaus]
						if get_lang == FORTRAN then
							pr @fact[3] === 0.0.to_var
						elsif get_lang == C
              code =<<EOF
							memset(&#{@fact[3]},0,sizeof(#{@fact[3]}));		
EOF
							get_output.print code
            end
						opn for_idime
 							pr @fact[3] === @fact[3] + @gpgve[@idime,@idime,@igaus]
						close for_idime
						pr @fact[3] === @fact[3] * @fact[1]
						opn for_inode
							opn for_jnode
								opn for_idime
									pr @idofv[1] === (@inode - 1) * @ndime + @idime
									pr @jdofv[1] === (@jnode - 1) * @ndime + @idime
									pr @elauu[@idofv[1],@jdofv[1]] === @elauu[@idofv[1],@jdofv[1]] + @fact[3] * @gpsha[@inode,@igaus] * @gpsha[@jnode,@igaus] 
								close for_idime
							close for_jnode
						close for_inode
					close for_igaus	
		}

		@block2 = lambda {
						opn for_igaus
							pr @fact[1] === 0.5.to_var * @gpden[@igaus] * @gpvol[@igaus]
							if get_lang == FORTRAN then
								pr @fact[3] === 0.0.to_var
							elsif get_lang == C
              	code =<<EOF
								memset(&#{@fact[3]},0,sizeof(#{@fact[3]})); 
EOF
								get_output.print code
            	end
							opn for_idime
								pr @fact[3] === @fact[3] + @gpgve[@idime,@idime,@igaus]
							close for_idime
							pr @fact[3] === @fact[3] * @fact[1]
							opn for_inode
								opn for_jnode
									opn for_idime
										pr @idofv[1] === (@inode - 1) * @ndime + @idime
                  	pr @jdofv[1] === (@jnode - 1) * @ndime + @idime
                  	pr @elauu[@idofv[1],@jdofv[1]] === @elauu[@idofv[1], @jdofv[1]] + @fact[3] * @gpsha[@inode,@igaus] * @gpsha[@jnode,@igaus]
									close for_idime
								close for_jnode
							close for_inode
						close for_igaus
						pr If(@kfl_regim_nsi == 3){
							opn for_itime
								 pr for2 = For(@igaus,1,@pgaus){
									if get_lang == FORTRAN then
                  	pr @gpveo === 0.0.to_var
                	elsif get_lang == C
                  	code =<<EOF
                  	memset(gpveo,0,sizeof(#{@gpveo.type.decl}) * #{@gpveo.dimension[0]});
EOF
										get_output.print code
									end
									opn for_inode
										opn for_idime
                    	pr @gpveo[@idime] ===  @gpveo[@idime] + @elvel[@idime,@inode,@itime] * @gpsha[@inode,@igaus]
                  	close for_idime
                	close for_inode
                	pr @fact[1] === 0.5.to_var * (@gpden[@igaus] - @densi[@igaus, @itime]) * @pabdf_nsi[@itime] * @dtinv_loc * @gpvol[@igaus]
                	opn for_inode
                  	opn for_idime
                    	pr @elrbu[@idime, @inode] === @elrbu[@idime,@inode] + @fact[1] * @gpsha[@inode, @igaus] * @gpveo[@idime]
                  	close for_idime
                	close for_inode
              	}
							close for_itime
						}
		}

		@block3 = lambda {
						opn for_igaus
							pr @fact[1] === @gpden[@igaus] * @gpvol[@igaus]
							if get_lang == FORTRAN then
								pr @fact[3] === 0.0.to_var
							elsif get_lang == C
                code =<<EOF
                memset(&#{@fact[3]},0,sizeof(#{@fact[3]})); 
EOF
                get_output.print code
              end
						  opn for_idime
	 							pr @fact[3] === @fact[3] + @gpgve[@idime,@idime,@igaus]
						  close for_idime
							pr @fact[3] === @fact[3] * @fact[1]
							opn for_inode
								opn for_jnode
									opn for_idime
										pr @idofv[1] === (@inode - 1) * @ndime + @idime
                  	pr @jdofv[1] === (@jnode - 1) * @ndime + @idime
                  	pr @elauu[@idofv[1],@jdofv[1]] === @elauu[@idofv[1], @jdofv[1]] + @fact[3] * @gpsha[@inode,@igaus] * @gpsha[@jnode,@igaus]
									close for_idime
								close for_jnode
							close for_inode 
							opn for_inode
								opn for_idime
									pr @idofv[1] === (@inode - 1)* @ndime + @idime
									opn for_jnode
										opn for_jdime
											pr @jdofv[1] === (@jnode - 1)* @ndime + @jdime
											pr @elauu[@idofv[1],@jdofv[1]] === @elauu[@idofv[1],@jdofv[1]] + @fact[1] * @gpsha[@inode,@igaus] * @gpcar[@idime, @jnode, @igaus] * @gpvel[@jdime, @igaus,1]
										close for_jdime
									close for_jnode
								close for_idime
							close for_inode
							opn for_inode
								opn for_idime
									opn for_jdime
										pr @elrbu[@idime,@inode] === @elrbu[@idime,@inode] + @fact[1] * @gpsha[@inode,@igaus] * @gpvel[@jdime,@igaus,1] * @gpgve[@idime,@jdime,@igaus] 
									close for_jdime
								close for_idime
							close for_inode
						close for_igaus
		}


    main_block = lambda{
			pr If(@cond1 => @block1, @cond2 => @block2, @cond3 => @block3)
    }
    construct(main_block)
  end
end

class Nest4 < Subroutine
    def initialize(options, functions = nil)
      super("nest4",options,functions)
      
      # ins
      @pgaus = Parameters.copy($pgaus,:in)
      @pnode = Parameters.copy($pnode,:in)
      @mnode = Parameters.copy($mnode,:in)
      @gpvis = Parameters.copy($gpvis,:in)
      @gpvol = Parameters.copy($gpvol,:in)
      @gpcar = Parameters.copy($gpcar,:in)
      @fvins_nsi = Parameters.copy($fvins_nsi,:in)
    	@ndime = Parameters.copy($ndime,:in)
    	@vsize = Parameters.copy($vsize,:in)

      # inouts
      @elauu = Parameters.copy($elauu,:inout)

      @args.push @vsize,@ndime,@pgaus,@pnode,@mnode,@gpvis,@gpvol,@gpcar,@fvins_nsi,@elauu
    end
  def generate
    # locals
    @igaus = Parameters.copy($igaus)
    @idime = Parameters.copy($idime)
    @jdime = Parameters.copy($jdime)
    @inode = Parameters.copy($inode)
    @jnode = Parameters.copy($jnode)
    @fact = Parameters.copy($fact)
    @idofv = Parameters.copy($idofv)
    @jdofv = Parameters.copy($jdofv)
    
    @locs.push @igaus, @idime, @jdime, @inode, @jnode, @fact, @idofv, @jdofv

    for_inode = For(@inode,1,@pnode)
    for_jnode = For(@jnode,1,@pnode)
    for_idime = For(@idime,1,@ndime)
    for_jdime = lambda{ |dim| return For(@jdime,1,dim) }
    
    exp_fact1 = @fact[2] === @gpvis[@igaus] * @gpvol[@igaus] * @gpcar[@idime,@jnode,@igaus]
    exp_idofv = @idofv[1] === (@inode-1)*@ndime + @idime
    
    block_elauu_1 = lambda{
      pr @jdofv[1] === (@jnode-1)*@ndime + @jdime
      pr @elauu[@idofv[1],@jdofv[1]] === @elauu[@idofv[1],@jdofv[1]] + @fact[2] * @gpcar[@jdime,@inode,@igaus]
    }
    
    exp_fact1_2 = @fact[2] === -2.0.to_var / 3.0.to_var * @gpvis[@igaus] * @gpvol[@igaus] * @gpcar[@idime,@inode,@igaus]
    
    block_elauu_2 = lambda{
      pr @jdofv[1] === (@jnode-1)*@ndime + @jdime
      pr @elauu[@idofv[1],@jdofv[1]] === @elauu[@idofv[1],@jdofv[1]] + @fact[2] * @gpcar[@jdime,@jnode,@igaus]
    }

    main_block = lambda{ 
 
      pr If(@fvins_nsi > 0.9){
        for1 = For(@igaus,1,@pgaus){ 
    
          opn for_inode
          opn for_idime
          pr exp_idofv
          opn for_jnode
          pr exp_fact1
          pr for_jdime.call(@ndime).unroll(@unroll), &block_elauu_1
          close for_jnode
          pr If(@fvins_nsi == 2.0){
            pr exp_fact1_2
            opn for_jnode
            pr for_jdime.call(@ndime).unroll(@unroll), &block_elauu_2
            close for_jnode
          }
          close for_idime
          close for_inode
        }
        
        pr for1
      }
    }              

    construct(main_block)
  end
end

class Nest5 < Subroutine
  def initialize(options, functions = nil)
    super("nest5",options,functions)

    # ins
    @pgaus = Parameters.copy($pgaus,:in)
    @pnode = Parameters.copy($pnode,:in)
    @gpvol = Parameters.copy($gpvol,:in)
    @gpden = Parameters.copy($gpden,:in)
    @dtinv_mod = Parameters.copy($dtinv_mod,:in)
    @gpsha = Parameters.copy($gpsha,:in)
    @elvel = Parameters.copy($elvel,:in)
    @kfl_lumped = Parameters.copy($kfl_lumped,:in)
    @ndime = Parameters.copy($ndime,:in)
    @vsize = Parameters.copy($vsize,:in)
    # inouts
    @gpveo = Parameters.copy($gpveo,:inout)
    @elrbu = Parameters.copy($elrbu,:inout)
    @elauu = Parameters.copy($elauu,:inout)
		# ADD TO ARGS
    @args.push @vsize,@ndime,@pgaus,@pnode,@gpvol,@gpden,@dtinv_mod,@gpsha,@elvel,@kfl_lumped,@gpveo,@elrbu,@elauu

    # locals 
    @idime = Parameters.copy($idime)
    @jdime = Parameters.copy($jdime)
    @igaus = Parameters.copy($igaus)
    @inode = Parameters.copy($inode) 
    @jnode = Parameters.copy($jnode)
    @fact = Parameters.copy($fact)
    @idof = Parameters.copy($idof)
    @jdof = Parameters.copy($jdof)

    @locs.push @idime, @jdime, @igaus, @inode, @jnode, @fact, @idof, @jdof
  end
  def generate


    exp_fact0 = @fact[1] === @gpvol[@igaus] * @gpden[@igaus] * @dtinv_mod            
    exp_fact1 = @fact[2] === @fact[1] * @gpsha[@inode,@igaus]
    exp_elrbu = @elrbu[@idime,@inode] === @elrbu[@idime,@inode] + @fact[2] * @elvel[@idime,@inode,2]
    for_elauu1 = For(@idime,1,@ndime){
      pr @idof[@idime] === (@inode-1)*@ndime+@idime
      pr @elauu[@idof[@idime],@idof[@idime]] === @elauu[@idof[@idime],@idof[@idime]] + @fact[2]
    }
    for_elrbu2 = For(@idime,1,@ndime){
      pr exp_elrbu
    }
    
    for_elauu2 = For(@jdime,1,3){
      pr @jdof[@jdime] === (@jnode-1)*@ndime+@jdime
      pr @elauu[@idof[@jdime],@jdof[@jdime]] === @elauu[@idof[@jdime],@jdof[@jdime]] - @fact[2] * @gpsha[@jnode,@igaus]
    }
    
    for_elrbu1 = For(@idime,1,3){
      pr @elrbu[@idime,@inode] ===  @elrbu[@idime,@inode] - @fact[2] * @gpveo[@idime]
      pr exp_elrbu
    }              
        block_kfl_lumped_1 = lambda{
          pr For(@igaus,1,@pgaus){
            if get_lang == FORTRAN then
              pr @gpveo === 0.0.to_var
            elsif get_lang == C
              code =<<EOF
              memset(gpveo,0,sizeof(#{@gpveo.type.decl}) * #{@gpveo.dimension[0]});
EOF
              get_output.print code
            end
            pr For(@inode,1,@pnode){
              pr For(@idime,1,3){
                pr @gpveo[@idime] === @gpveo[@idime] + @elvel[@idime,@inode,2] * @gpsha[@inode,@igaus] 
              }.unroll(@unroll)

            }
            pr exp_fact0
            pr For(@inode,1,@pnode){
              pr exp_fact1
              if @unroll then
                pr for_elauu1.unroll
                pr for_elrbu1.unroll
              else
                pr for_elauu1
                pr for_elrbu1
              end
              pr For(@jnode,1,@pnode){
                pr for_elauu2.unroll(@unroll)
              }
            }
          }
        }
    block_kfl_lumped_2 = lambda{
        pr For(@igaus,1,@pgaus){
          pr exp_fact0
          pr For(@inode,1,@pnode){
            pr exp_fact1
            pr for_elauu1.unroll(@unroll)
            pr for_elrbu2.unroll(@unroll)
          }
        }
      }
    if_kfl_lumped_1 = @ndime == 2 ? lambda{puts "PREGUNTAR A MATIAS QUE LO PROGRAME"} : block_kfl_lumped_1
    if_kfl_lumped_2 = block_kfl_lumped_2
    
    main_block = lambda{ 
      pr If(@kfl_lumped == 1 => if_kfl_lumped_1, @kfl_lumped == 2 => if_kfl_lumped_2) 
    }      
    construct(main_block)
  end
end

class Nest6 < Subroutine
  def initialize(options, functions = nil)
    super("nest6",options,functions)
    
    # ins
    @mnode = Parameters.copy($mnode,:in)
    @pgaus = Parameters.copy($pgaus,:in)
    @pnode = Parameters.copy($pnode,:in)
    @gpvol = Parameters.copy($gpvol,:in)
    @gpsha = Parameters.copy($gpsha,:in)
    @gpcar = Parameters.copy($gpcar,:in)
    @ndime = Parameters.copy($ndime,:in)
    @vsize = Parameters.copy($vsize,:in)
    # inouts
    @elapu = Parameters.copy($elapu,:inout)
    @elaup = Parameters.copy($elaup,:inout)
		# ADD TO ARGS
    @args.push @vsize,@ndime,@mnode,@pgaus,@pnode,@gpvol,@gpsha,@gpcar,@elapu,@elaup

    # locals  
    @igaus = Parameters.copy($igaus)
    @inode = Parameters.copy($inode)
    @jnode = Parameters.copy($jnode)
    @idime = Parameters.copy($idime)
    @jdime = Parameters.copy($jdime)
    @idof = Parameters.copy($idof)
    @fact = Parameters.copy($fact)
		# ADD TO LOCS
    @locs.push @igaus, @inode, @jnode, @idime, @jdime, @idof, @fact 
  end

  def generate
    for1 = For(@igaus,1,@pgaus){
      pr For(@inode,1,@pnode){
        pr For(@idime,1,@ndime){
          pr @idof[@idime] === (@inode-1)*@ndime + @idime
        }.unroll(@unroll)
    
        pr For(@jnode,1,@pnode){
          pr @fact[1] === @gpvol[@igaus] * @gpsha[@jnode,@igaus]
          pr For(@jdime,1,@ndime){
            pr @fact[@jdime+1] === @fact[1] * @gpcar[@jdime,@inode,@igaus]
            pr @elapu[@jnode,@idof[@jdime]] === @elapu[@jnode,@idof[@jdime]] + @fact[@jdime+1]
          }.unroll(@unroll)
    
          pr For(@jdime,1,@ndime){
            pr @elaup[@idof[@jdime],@jnode] === @elaup[@idof[@jdime],@jnode] - @fact[@jdime+1]
          }.unroll(@unroll)
        }
      }
    }
    
    main_block = lambda{
      pr for1
    }
    construct(main_block)
  end
end

# NEST MODIFIED
class Nest7 < Subroutine
  def initialize(options, functions = nil)
    super("nest7",options,functions)
    
    # ins
    @kfl_stabi_nsi = Parameters.copy($kfl_stabi_nsi,:in)
    @pgaus = Parameters.copy($pgaus, :in)
    @pnode = Parameters.copy($pnode, :in)
    @gpsp1_p = Parameters.copy($gpsp1_p, :in)
    @wgrgr = Parameters.copy($wgrgr, :in)
    @gpvol = Parameters.copy($gpvol, :in)
    @nsi_galerkin = Parameters.copy($nsi_galerkin, :in)
    @nsi_algebraic_split_oss = Parameters.copy($nsi_algebraic_split_oss, :in)
    @vsize = Parameters.copy($vsize,:in)
    # inouts
    @elapp = Parameters.copy($elapp, :inout)
		# ADD TO ARGS
    @args.push @vsize,@kfl_stabi_nsi,@pgaus,@pnode,@gpsp1_p,@wgrgr,@gpvol,@nsi_galerkin,@nsi_algebraic_split_oss,@elapp

    # locals
    @igaus = Parameters.copy($igaus)
    @inode = Parameters.copy($inode)
    @jnode = Parameters.copy($jnode)
    @fact = Parameters.copy($fact)
		# ADD TO LOCS
    @locs.push @igaus, @inode, @jnode, @fact
  end

  def generate
    main_block = lambda{
      pr If( And.basic_usage(@kfl_stabi_nsi != @nsi_galerkin , @kfl_stabi_nsi != @nsi_algebraic_split_oss) ){
        pr For(@igaus,1,@pgaus){
          pr For(@inode,1,@pnode){
            pr For(@jnode,@inode+1,@pnode){
              pr @fact[2]             === @gpsp1_p[@igaus] * @wgrgr[@jnode,@inode,@igaus] * @gpvol[@igaus]
              pr @elapp[@jnode,@inode] === @elapp[@jnode,@inode] + @fact[2]
              pr @elapp[@inode,@jnode] === @elapp[@inode,@jnode] + @fact[2]
            }
            pr @fact[2]             === @gpsp1_p[@igaus] * @wgrgr[@inode,@inode,@igaus] * @gpvol[@igaus]
            pr @elapp[@inode,@inode] === @elapp[@inode,@inode] + @fact[2]
          }
        }
      } 
    }      

    construct(main_block)
  end
end

class Nest8 < Subroutine
  def initialize(options, functions = nil)
    super("nest8",options,functions)

    # ins
    @pnode = Parameters.copy($pnode, :in)
    @pgaus = Parameters.copy($pgaus, :in)
    @penal_nsi = Parameters.copy($penal_nsi, :in)
    @gpvol = Parameters.copy($gpvol, :in)
    @gpsha = Parameters.copy($gpsha, :in)
    @elpre = Parameters.copy($elpre, :in)
    @vsize = Parameters.copy($vsize,:in)
    # inouts
    @elapp = Parameters.copy($elapp, :inout)
    @elrbp = Parameters.copy($elrbp, :inout)
		# ADD TO ARGS
    @args.push @vsize,@pnode,@pgaus,@penal_nsi,@gpvol,@gpsha,@elpre,@elapp,@elrbp

    # locals
    @fact = Parameters.copy($fact)
    @igaus = Parameters.copy($igaus)
    @inode = Parameters.copy($inode)
		# ADD TO LOCS
    @locs.push @fact, @igaus, @inode
  end  

  def generate
    main_block = lambda {
      pr For(@igaus,1,@pgaus){
        pr @fact[2] === @penal_nsi * @gpvol[@igaus]
        pr For(@inode,1,@pnode){
          pr @elapp[@inode,@inode] === @elapp[@inode,@inode] + @fact[2] * @gpsha[@inode,@igaus]
          pr @elrbp[@inode]       === @elrbp[@inode]       + @fact[2] * @gpsha[@inode,@igaus] * @elpre[@inode,1]
        }
      }
    }
    construct(main_block)
  end
end

class Nest9 < Subroutine
  def initialize(options, functions = nil)
    super("nest9",options,functions)
    
    # ins
    @pgaus = Parameters.copy($pgaus, :in)
    @pnode = Parameters.copy($pnode, :in)
    @elvel = Parameters.copy($elvel, :in)
    @agrau = Parameters.copy($agrau, :in)
    @gpsp1 = Parameters.copy($gpsp1, :in)
    @kfl_limit_nsi = Parameters.copy($kfl_limit_nsi, :in)
    @kfl_stabi_nsi = Parameters.copy($kfl_stabi_nsi, :in)
    @nsi_galerkin = Parameters.copy($nsi_galerkin, :in)
    @nsi_algebraic_split_oss = Parameters.copy($nsi_algebraic_split_oss, :in)
    @ndime = Parameters.copy($ndime,:in)
    @vsize = Parameters.copy($vsize,:in)
    # inouts
    @gpvep = Parameters.copy($gpvep, :inout)
		# ADD TO ARGS
    @args.push @vsize,@ndime,@pgaus,@pnode,@elvel,@agrau,@gpsp1,@kfl_limit_nsi,@kfl_stabi_nsi,@nsi_galerkin,@nsi_algebraic_split_oss,@gpvep

    # local
    @igaus = Parameters.copy($igaus)
    @idime = Parameters.copy($idime)
    @inode = Parameters.copy($inode)
    @c1 = Parameters.copy($c1)
    @c2 = Parameters.copy($c2)
    @c3 = Parameters.copy($c3)
    @c4 = Parameters.copy($c4)
    @beta = Parameters.copy($beta)
    @alpha = Parameters.copy($alpha)
		# ADD TO LOCS
    @locs.push @igaus,@idime,@inode,@c1,@c2,@c3,@c4,@beta,@alpha

  end
  def generate

    register_funccall("epsilon") if get_lang == FORTRAN
    
    tanh = lambda{|x|
      if get_lang == FORTRAN or @vector_length < 2 then
        return Tanh(x)
      else
        return @functions[:tanh].call(x)
      end
    }
    
    min = lambda{|x,y|
      if get_lang == FORTRAN then 
        return Min( x, y )
      else
        if @vector_length > 1 then
          return @functions[:min].call(x,y)
        else
          return Ternary(x < y, x, y)
        end
      end
    }
    for_igaus = For(@igaus,1,@pgaus){
     if get_lang == FORTRAN then
      pr @c1 === 0.0.to_var 
      pr @c2 === 0.0.to_var 
      pr @c3 === 0.0.to_var 
     elsif get_lang == C
      code =<<EOF
			memset(&#{@c1},0,sizeof(#{@c1}));
			memset(&#{@c2},0,sizeof(#{@c2}));
			memset(&#{@c3},0,sizeof(#{@c3}));
EOF
      get_output.print code
     end
     pr For(@idime,1,@ndime){
     if get_lang == FORTRAN then
       pr @c4 === 0.0.to_var
     elsif get_lang == C
       code =<<EOF
			 memset(&#{@c4},0,sizeof(#{@c4}));
EOF
       get_output.print code
     end
     pr For(@inode,1,@pnode){
        pr @c4 === @c4 + @agrau[@inode,@igaus] * @elvel[@idime,@inode,1]
     }
     pr @c4 === @gpsp1[@igaus] * @c4
     # Exponential intrinsics only works with ICC
     pr @c1 === @c1 + ( @gpvep[@idime,@igaus] - @c4 ) * ( @gpvep[@idime,@igaus] - @c4 )
     pr @c3 === @c3 + @gpvep[@idime,@igaus] * @gpvep[@idime,@igaus]
     pr @c2 === @c2 + @c4 * @c4
    }
    
    pr @c3   === Sqrt( @c2 ) + Sqrt( @c3 )  
    pr @c1   === Sqrt( @c1 )
    if get_lang == FORTRAN then
      pr @beta === @c1 / ( @c3 + epsilon(1.0.to_var) )
    elsif get_lang == C
      pr @beta === @c1 / ( @c3 + Real(:DBL_EPSILON) );
    end
    pr If(@kfl_limit_nsi == 1 => lambda{
      pr @alpha === min.call( Set(1.0.to_var, @alpha), 2.0.to_var * ( 1.0.to_var  - @beta ) )
    }, @kfl_limit_nsi == 2 => lambda{
      pr @alpha === 0.5 * ( tanh.call( 20.0 * ( @beta - 0.8 ) ) + 1.0 )
    })              
    pr For(@idime,1,@ndime){
      pr @gpvep[@idime,@igaus] === @alpha * @gpvep[@idime,@igaus]
    }
   }

   if get_lang == C then
   code =<<EOF
      memset(gpvep,0,sizeof(#{@gpvep.type.decl}) * #{@ndime} * pgaus);
EOF
   end

   main_block = lambda{
     pr If( Or.basic_usage( Or.basic_usage(@kfl_limit_nsi == -1, @kfl_stabi_nsi == @nsi_galerkin), @kfl_stabi_nsi == @nsi_algebraic_split_oss ) => lambda{
     if get_lang == FORTRAN then
         pr @gpvep === 0.0
     elsif get_lang == C
         get_output.print code
     end           
     }, @kfl_limit_nsi > 0 => lambda{ pr for_igaus } )
   }
   construct(main_block)
  end
end

class Nest10 < Subroutine
  def initialize(options, functions = nil)
    super("nest10",options,functions)
    
    # ins
    @pgaus = Parameters.copy($pgaus, :in)
    @gpsp1_p = Parameters.copy($gpsp1_p, :in)
    @gprhs = Parameters.copy($gprhs, :in)
    @gpden = Parameters.copy($gpden, :in)
    @dtsgs = Parameters.copy($dtsgs, :in)
    @gpsp1_v = Parameters.copy($gpsp1_v, :in)
    @gpsgs = Parameters.copy($gpsgs, :in)
    @kfl_sgsti_nsi = Parameters.copy($kfl_sgsti_nsi, :in)
    @kfl_stabi_nsi = Parameters.copy($kfl_stabi_nsi, :in)
    @nsi_galerkin = Parameters.copy($nsi_galerkin, :in)
    @nsi_algebraic_split_oss = Parameters.copy($nsi_algebraic_split_oss, :in)
    @ndime = Parameters.copy($ndime,:in)
    @vsize = Parameters.copy($vsize,:in)
    # inouts
    @gpgrp = Parameters.copy($gpgrp, :inout)
    @gpvep = Parameters.copy($gpvep, :inout)
		# ADD TO ARGS
    @args.push @vsize,@ndime,@pgaus,@gpsp1_p,@gprhs,@gpden,@dtsgs,@gpsp1_v,@gpsgs,@kfl_sgsti_nsi,@kfl_stabi_nsi,@nsi_galerkin,@nsi_algebraic_split_oss,@gpgrp,@gpvep

    # local
    @igaus = Parameters.copy($igaus)
    @idime = Parameters.copy($idime)
    @fact = Parameters.copy($fact)
    @fact1_p = Parameters.copy($fact1_p)
		# ADD TO LOCS
    @locs.push @igaus,@idime,@fact,@fact1_p

  end
  def generate
		for_igaus = []
    for_igaus[0] = For(@igaus,1,@pgaus){
        pr For(@idime,1,@ndime){
            pr @gpgrp[@idime,@igaus] === @gpgrp[@idime,@igaus] + @gpsp1_p[@igaus] * @gprhs[@idime,@igaus]
        }.unroll(@unroll)
    }
              
    for_igaus[1] = For(@igaus,1,@pgaus){
        pr @fact[2]     === @gpden[@igaus] * @dtsgs * @gpsp1_v[@igaus]
        pr @fact1_p[1]  === @gpden[@igaus] * @dtsgs * @gpsp1_p[@igaus]
        pr For(@idime,1,@ndime){
             pr @gpvep[@idime,@igaus] === @gpvep[@idime,@igaus] + @fact[2]    * @gpsgs[@idime,@igaus,2]
             pr @gpgrp[@idime,@igaus] === @gpgrp[@idime,@igaus] + @fact1_p[1] * @gpsgs[@idime,@igaus,2]
        }.unroll(@unroll)
    }

    block = lambda{
        pr for_igaus[0]
        pr If(@kfl_sgsti_nsi == 1){
             pr for_igaus[1]
        }
    }

    main_block = lambda{
        pr If( Or.basic_usage(@kfl_stabi_nsi == @nsi_galerkin, @kfl_stabi_nsi == @nsi_algebraic_split_oss ) => lambda{ 
             if get_lang == FORTRAN then
                   pr @gpgrp === 0.0 
             elsif get_lang == C
                   code =<<EOF
                   memset(gpgrp,0,sizeof(#{@gpgrp.type.decl}) * #{@ndime} * pgaus);
EOF
             end
             get_output.print code
        }, :else => block)
    }
    construct(main_block)
  end
end

class Nest11 < Subroutine
  def initialize(options, functions = nil)
    super("nest11",options,functions)
    
    # ins
    @mnode = Parameters.copy($mnode, :in)
    @pgaus = Parameters.copy($pgaus, :in)
    @gpden = Parameters.copy($gpden, :in)
    @dtinv_mod = Parameters.copy($dtinv_mod, :in)
    @nbdfp_nsi = Parameters.copy($nbdfp_nsi, :in)
    @gpvel = Parameters.copy($gpvel, :in)
    @pnode = Parameters.copy($pnode, :in)
    @gpvol = Parameters.copy($gpvol, :in)
    @gpsha = Parameters.copy($gpsha, :in)
    @agrau = Parameters.copy($agrau, :in)
    @gprhc = Parameters.copy($gprhc, :in)
    @gpcar = Parameters.copy($gpcar, :in)
    @gpvep = Parameters.copy($gpvep, :in)
    @gpgrp = Parameters.copy($gpgrp, :in)
    @ndime = Parameters.copy($ndime,:in)
    @vsize = Parameters.copy($vsize,:in)
    # inouts
    @gprhs = Parameters.copy($gprhs, :inout)
    @elrbu = Parameters.copy($elrbu, :inout)
    @elrbp = Parameters.copy($elrbp, :inout)
    @pabdf_nsi = Parameters.copy($pabdf_nsi, :inout)
		# ADD TO ARGS
    @args.push @vsize,@ndime,@mnode,@pgaus,@gpden,@dtinv_mod,@nbdfp_nsi,@gpvel,@pnode,@gpvol,@gpsha,@agrau,@gprhc,@gpcar,@gpvep,@gpgrp,@gprhs,@elrbu,@elrbp,@pabdf_nsi

    # locals
    @igaus = Parameters.copy($igaus)
    @fact = Parameters.copy($fact)
    @idime = Parameters.copy($idime)
    @itime = Parameters.copy($itime)
    @inode = Parameters.copy($inode)
		# ADD TO LOCS
    @locs.push @igaus,@fact,@idime,@itime,@inode

  end

  def generate

      for_igaus = For(@igaus,1,@pgaus){
        pr @fact[5] === @gpden[@igaus] * @dtinv_mod
        pr For(@itime,2,@nbdfp_nsi){
          pr For(@idime,1,@ndime){
            #pr @gprhs[@idime,@igaus] === @gprhs[@idime,@igaus] - @functions[:pabdf_nsi].call(@itime) * @fact[5] * @gpvel[@idime,@igaus,@itime]
            pr @gprhs[@idime,@igaus] === @gprhs[@idime,@igaus] - @pabdf_nsi[@itime] * @fact[5] * @gpvel[@idime,@igaus,@itime]
          }.unroll(@unroll)
        }
        pr For(@inode,1,@pnode){
          pr @fact[2] === @gpvol[@igaus] * @gpsha[@inode,@igaus]  # ( f + rho*u^n/dt , v )
          pr @fact[4] === @gpvol[@igaus] * @agrau[@inode,@igaus]  # ( rho * a.grad(v) , P1' ) 
          pr For(@idime,1,@ndime){
            pr @elrbu[@idime,@inode] === @elrbu[@idime,@inode] + @fact[2] * @gprhs[@idime,@igaus] + @fact[4] * @gpvep[@idime,@igaus]
          }.unroll(@unroll)
          pr @elrbp[@inode] === @elrbp[@inode] + @gpvol[@igaus] * @gpsha[@inode,@igaus] * @gprhc[@igaus]  # ( rhs, q )
          pr For(@idime,1,@ndime){
            pr @elrbp[@inode] === @elrbp[@inode] + @gpvol[@igaus] * @gpcar[@idime,@inode,@igaus] * @gpgrp[@idime,@igaus] # ( P2' , grad(q) ) 
          }.unroll(@unroll)
        }
      }
    
    main_block = lambda{
      pr for_igaus
    }
    construct(main_block)
    end
end

class Nest12 < Subroutine
            def initialize(options, functions = nil)
              super("nest12",options,functions)
 
              # ins
              @mnode = Parameters.copy($mnode, :in)
              @pgaus = Parameters.copy($pgaus, :in)
              @pnode = Parameters.copy($pnode, :in)
              @gpcar = Parameters.copy($gpcar, :in)
              @gpvol = Parameters.copy($gpvol, :in)
              @gpsha = Parameters.copy($gpsha, :in)
              @gpcar_bub = Parameters.copy($gpcar_bub, :in)
              @pbubl = Parameters.copy($pbubl, :in)
              @kfl_stabi_nsi = Parameters.copy($kfl_stabi_nsi, :in)
              @kfl_press_nsi = Parameters.copy($kfl_press_nsi, :in)
              @penal_nsi = Parameters.copy($penal_nsi, :in)
              @elbub = Parameters.copy($elbub, :in)
              @gprhc = Parameters.copy($gprhc, :in)
              @gpsha_bub = Parameters.copy($gpsha_bub, :in)
              @nsi_galerkin = Parameters.copy($nsi_galerkin, :in)
              @nsi_algebraic_split_oss = Parameters.copy($nsi_algebraic_split_oss, :in)
    					@ndime = Parameters.copy($ndime,:in)
    					@vsize = Parameters.copy($vsize,:in)
              # inouts
              @elauq = Parameters.copy($elauq, :inout)
              @elaqu = Parameters.copy($elaqu, :inout)
              @elapq = Parameters.copy($elapq, :inout)
              @elaqp = Parameters.copy($elaqp, :inout)
              @elaqq = Parameters.copy($elaqq, :inout)
              @elrbq = Parameters.copy($elrbq, :inout)
							# ADD TO ARGS
              @args.push @vsize,@ndime,@mnode,@pgaus,@pnode,@gpcar,@gpvol,@gpsha,@gpcar_bub,@pbubl,@kfl_stabi_nsi,@kfl_press_nsi,@penal_nsi,@elbub,@gprhc,@gpsha_bub,@elauq,@elaqu,@elapq,@elaqp,@elaqq,@elrbq, @nsi_galerkin,@nsi_algebraic_split_oss

              # locals
              @fact = Parameters.copy($fact)
              @igaus  = Parameters.copy($igaus)
              @idof = Parameters.copy($idof)
              @idime = Parameters.copy($idime)
              @inode = Parameters.copy($inode)
							# ADD TO LOCS
	      			@locs.push @fact,@igaus,@idof,@idime,@inode

            end

            def generate

              register_funccall("stop") if get_lang == FORTRAN
              register_funccall("maxval") if get_lang == FORTRAN

              exp_elauq = []
              exp_elauq[0] = "@elauq[@idof[@ndime],1] === @elauq[@idof[@ndime],1] - @fact[2] * @gpcar[@idime,@inode,@igaus]"
              exp_elauq[1] = "@elauq[@idof[@ndime],1] === @elauq[@idof[@ndime],1] + @gpvol[@igaus] * @gpsha[@inode,@igaus] * @gpcar_bub[@idime,@igaus]"
              
              block_igaus = lambda{|exp|
                return For(@igaus,1,@pgaus){
                  pr @fact[2] === @gpvol[@igaus] * @gpsha_bub[@igaus]
                  pr For(@inode,1,@pnode){
                    pr For(@idime,1,@ndime){
                      pr @idof[@ndime] === (@inode-1)*@ndime + @idime
                      pr eval exp
                      pr @elaqu[1,@idof[@ndime]] === @elaqu[1,@idof[@ndime]] + @fact[2] * @gpcar[@idime,@inode,@igaus]
                    }.unroll(@unroll)
                  }
                }
              }
              
              if_kfl_press_nsi_1 = block_igaus.call(exp_elauq[0])
              if_kfl_press_nsi_else = block_igaus.call(exp_elauq[1])

              call_maxval = lambda{|x|
                if get_lang == FORTRAN
                  return maxval(x)
                else
                  return @functions[:maxval].call(x)
                end
              }

              main_block = lambda{
                pr If(call_maxval.call(@pbubl) == 1){
                  pr If( Or.basic_usage( @kfl_stabi_nsi != @nsi_galerkin, @kfl_stabi_nsi == @nsi_algebraic_split_oss)){
                    if get_lang == FORTRAN then
                      puts "BUBBLE NOT CODED FOR SPLIT OSS"
                    end
                  }
                  
                  if get_lang == FORTRAN then
                    pr @elauq === 0.0
                    pr @elapq === 0.0
                    pr @elaqu === 0.0
                    pr @elaqp === 0.0
                    pr @elaqq === 0.0
                    pr @elrbq === 0.0
                  elsif get_lang == C
                    code =<<EOF
                    memset(elauq, 0, sizeof(#{@elauq.type.decl}) * pnode * #{@ndime});
                    memset(elapq, 0, sizeof(#{@elapq.type.decl}) * pnode);
                    memset(elaqu, 0, sizeof(#{@elaqu.type.decl}) * pnode * #{@ndime});
                    memset(elaqp, 0, sizeof(#{@elaqp.type.decl}) * pnode);
                    memset(elaqq, 0, sizeof(#{@elaqq.type.decl}));
                    memset(elrbq, 0, sizeof(#{@elrbq.type.decl}));
EOF
                    get_output.print code
                  end
                  pr If(@kfl_press_nsi == 1 => lambda{ pr if_kfl_press_nsi_1}, :else => lambda{ pr if_kfl_press_nsi_else })
                  
                  pr For(@igaus,1,@pgaus){
                    pr @elaqq[1,1] === @elaqq[1,1] + @gpvol[@igaus] * @gpsha_bub[@igaus] * @penal_nsi
                    pr @elrbq[1]   === @elrbq[1]   + @gpvol[@igaus] * @gpsha_bub[@igaus] * @penal_nsi * @elbub
                    pr @elrbq[1]   === @elrbq[1]   + @gpvol[@igaus] * @gpsha_bub[@igaus] * @gprhc[@igaus]
                  }
                }
              }
              
              construct(main_block)
            end
          end



#.. CLASS KSplitBoast ..#

class KSplitBoast < KSplitOss

  def generate

    #.. FUNCTIONS ..#
    register_funccall("maxval") if get_lang == FORTRAN
    register_funccall("epsilon") if get_lang == FORTRAN
    register_funccall("stop") if get_lang == FORTRAN

    # maxval
    unless get_lang == FORTRAN then
      #x2 = Int("x", :dir => :in, :vector_length => @opts[:vector_length], :dim => [Dim(1)])
      x2 = Int("x", :dir => :in, :vector_length => @opts[:vector_length])
      y2 = Int("y")
      @p_maxval = Procedure("maxval",[x2], :return => y2){
        if @opts[:vector_length] > 1 then
          decl i = Int("i")
          decl a = Int("a", :dim => [Dim(@opts[:vector_length])], :allocate => true)
    
          #pr a[1] === x2.dereference
          pr a[1] === x2
          pr y2 === a[1]
    
          pr For(i,2,@opts[:vector_length]){
            pr If(a[i] > y2){
              pr y2 === a[i]
            }
          }
        else
          pr y2 === x2
        end
      }            
    end     
    # Tanh     
    unless get_lang == FORTRAN then
      x3 = Real("x", :dir => :in, :vector_length => @opts[:vector_length])
      y3 = Real("y", :vector_length => @opts[:vector_length])
    
      @p_tanh = Procedure("Tanh", [x3], :return => y3){
        decl i = Int("i")
        decl tmp1 = Real("tmp1", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        decl tmp2 = Real("tmp2", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        
        pr tmp1[1] === x3
    
        pr For(i,1,@opts[:vector_length]){
          pr tmp2[i] === Tanh(tmp1[i])
        }
        pr y3 === tmp2[1] 
      }
    end
		# Sqrt
    unless get_lang == FORTRAN then
      x4 = Real("x", :dir => :in, :vector_length => @opts[:vector_length])
      y4 = Real("y", :vector_length => @opts[:vector_length])
    
      @p_sqrt = Procedure("Sqrt", [x4], :return => y4){
        decl i = Int("i")
        decl tmp1 = Real("tmp1", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        decl tmp2 = Real("tmp2", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        
        pr tmp1[1] === x4
    
        pr For(i,1,@opts[:vector_length]){
          pr tmp2[i] === Sqrt(tmp1[i])
        }
        pr y4 === tmp2[1] 
      }
    end
    # Min
    unless get_lang == FORTRAN then
      x5 = Real("x", :dir => :in, :vector_length => @opts[:vector_length])
      y5 = Real("y", :dir => :in, :vector_length => @opts[:vector_length])
      z5 = Real("z", :vector_length => @opts[:vector_length])
    
      @p_min = Procedure("Min", [x5,y5], :return => z5){
        decl i = Int("i")
        decl tmp1 = Real("tmp1", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        decl tmp2 = Real("tmp2", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        decl tmp3 = Real("tmp3", :dim => [Dim(@opts[:vector_length])], :allocate => true)
        
        pr tmp1[1] === x5
        pr tmp2[1] === y5
    
        pr For(i,1,@opts[:vector_length]){
          pr tmp3[i] === Ternary(tmp1[i] < tmp2[i], tmp1[i], tmp2[i])
        }
        pr z5 === tmp3[1] 
      }
    end              


		#.. NESTS ..#
    #@opts[:ndime] = @ndime
    nests = {}
    @opts[:nests].each{|i|
      case i
      when 9
        functions = (get_lang == FORTRAN ? nil :{:tanh => @p_tanh, :min => @p_min})
      when 12
        functions = (get_lang == FORTRAN ? nil : {:maxval => @p_maxval})
      else
        functions = nil
      end
      eval "nests[:nest#{i}] = Nest#{i}::new(@opts,functions)"
      eval "nests[:nest#{i}].generate"
    }
  
    includes = ["immintrin.h"]
    includes.push("string.h", "math.h", "float.h") if get_lang == C
    @kernel = CKernel::new( :includes => includes )
    functions = []
    functions.push @p_maxval unless get_lang == FORTRAN
    @kernel.procedure = declare_procedure("nsi_element_assembly_split_oss_boast",functions.flatten)


    # Print subroutines           
    unless get_lang == FORTRAN then
      pr @p_maxval
      pr @p_min
      pr @p_tanh
    end
 
		# Print all nests
    nests.values.each{|n|
			pr n.code 
		} unless @opts[:inline] == :included




    # Main procedure body
    opn @kernel.procedure


 
if get_lang == FORTRAN
	pr If($nsi_fractional_step_int == 1 => lambda{
  	pr @dtinv_mod === 0.0
	}, :else => lambda{
  	pr @dtinv_mod === @dtinv_loc
	})

  pr @gpsp1_p === @gpsp1
  pr @gpsp1_v === @gpsp1
  pr @gpsp2_v === @gpsp2

  pr If( @kfl_nota1_nsi == 1 => lambda{ pr @gpsp1_v === 0.0 } )

	pr If( @kfl_stabi_nsi == $nsi_algebraic_split_oss => lambda{ 
		pr @gpsp1_p === 0.0
    pr @gpsp1_v === 0.0
    pr @gpsp2_v === 0.0
	})

  pr @elrbp === 0.0
  pr @elrbu === 0.0
  pr @elapp === 0.0
  pr @elauu === 0.0
  pr @elaup === 0.0
  pr @elapu === 0.0
  pr @agrau === 0.0
  pr @wgrgr === 0.0 
elsif get_lang == C
  code =<<EOF
	if(NSI_FRACTIONAL_STEP_int == 1)
		memset(&#{@dtinv_mod},0,sizeof(#{@dtinv_mod})); 
	else
		memcpy(&#{@dtinv_mod},&#{@dtinv_loc},sizeof(#{@dtinv_mod})); 

  memcpy(gpsp1_p, gpsp1, sizeof(gpsp1_p));
  memcpy(gpsp1_v, gpsp1, sizeof(gpsp1_v));
  memcpy(gpsp2_v, gpsp2, sizeof(gpsp2_v));

  if (kfl_nota1_nsi == 1) memset(gpsp1_v, 0, sizeof(gpsp1_v));

  if(kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS){
  	memset(gpsp1_p, 0, sizeof(#{@gpsp1_p.type.decl}) * pgaus);
  	memset(gpsp1_v, 0, sizeof(#{@gpsp1_v.type.decl}) * pgaus);
  	memset(gpsp2_v, 0, sizeof(#{@gpsp2_v.type.decl}) * pgaus);
	}

  memset(elrbp, 0, sizeof(#{@elrbp.type.decl}) * pnode);
  memset(elrbu, 0, sizeof(#{@elrbu.type.decl}) * #{@ndime} * pnode);
  memset(elapp, 0, sizeof(#{@elapp.type.decl}) * pnode * pnode);
  memset(elauu, 0, sizeof(#{@elauu.type.decl}) * pnode * #{@ndime} * pnode * #{@ndime});
  memset(elaup, 0, sizeof(#{@elaup.type.decl}) * pnode * #{@ndime} * pnode);
  memset(elapu, 0, sizeof(#{@elapu.type.decl}) * pnode * pnode * #{@ndime});
  memset(agrau, 0, sizeof(#{@agrau.type.decl}) * pnode * pgaus);
  memset(wgrgr, 0, sizeof(#{@wgrgr.type.decl}) * pnode * pnode * pgaus);
EOF
  get_output.print code
end

      # Either call subroutine or paste their code
      nests.values.each{|n| n.call } 

    close @kernel.procedure
  end        
end
