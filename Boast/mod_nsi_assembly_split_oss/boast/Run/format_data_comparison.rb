require 'yaml'
require 'pp'
require 'csv'


# USAGE : ruby format_data_comparison.rb ../experiments/../Data.yaml test.cvs
input = ARGV[0]
output = ARGV[1]

#header = [:kernel,:vector_length,:inline,:unroll,:FCFLAGS,:time]
header = []
body = []

h = YAML::load(File::open(input).read)
h.each{|k1,v1|
		e = v1[:time].min
    t = []
		t.push k1[:kernel]
		t.push k1[:vector_length]
		t.push k1[:inline]
		t.push k1[:unroll]
	  t.push k1[:FCFLAGS]
    t.push(e)
    body.push(t)
}

CSV.open(output, "w"){ |f|
  f << header
  body.each{ |e|
    f << e
  }
}

