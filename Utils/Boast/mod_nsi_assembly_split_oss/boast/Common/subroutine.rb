require_relative './Parameters.rb'

class Subroutine

  include Parameters
  attr_accessor :code
  def initialize(name,options,functions)
    opts = {:vector_length => 1,:unroll => false, :inline => :included}
    opts.update(options)
    @name = name
    @functions = functions
    @args = []
    @locs = []
    @code = nil
    @vector_length = options[:vector_length]
    @usage = opts[:inline]
    @unroll = opts[:unroll]
    #@ndime = opts[:ndime]
    @nb_original_param = instance_variables.length + 1
  end

  def construct(block)
    if @usage == :included then
      @code = block
    else
      functions = @functions.nil? ? nil : @functions.values
      inlined = @usage == :inlined ? true : false
      @code = Procedure(@name, @args, :inline => inlined, :functions => functions, :locals => @locs , &block)
    end
  end

  def call
    if @usage == :included then
      @code.call
    else
      pr @code.call( *@args )
    end
  end

end
