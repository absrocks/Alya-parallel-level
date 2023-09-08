# Running Alya using BOAST (MN4) {#boast-general}

This page explains how to run Alya by taking advantage of kernel optimization provided by BOAST.
[BOAST](https://github.com/Nanosim-LIG/boast) is a framework to metaprogram, benchmark and validate computing kernels. It enables a "semi-automatic" optimization of C and Fortran kernels, depending on the platform and compiler options that are used.
This tutorial is aimed at Marenotrum IV, but you can easily adapt it to another platform. You need the last version of Alya's sources, as well as ruby installed, and the following ruby modules:

    'narray', '~> 0.6.0', '>=0.6.0.8'
    'narray_ffi', '~> 1.2', '>=1.2.0'
    'opencl_ruby_ffi', '~> 1.3', '>=1.3.2'
    'systemu', '~> 2', '>=2.2.0'
    'os', '~> 0.9', '>=0.9.6'
    'PAPI', '~> 1.0', '>=1.0.0'
    'hwloc', '~> 0.3', '>=0.3.0'
    'ffi', '~> 1.9', '>=1.9.3'
    'rgl', '~> 0.5', '>=0.5.1'

## Install BOAST on MN4

In order to take advantage of BOAST, you need to install on your MN4 account.
First, you have to retrieved BOAST's git repository.

    $ git clone https://github.com/Nanosim-LIG/boast.git

And then copy it to your MN account:

    $ scp -r bscxxxxx@mn1.bsc.es:~

Connect to MN:

    $ ssh bscxxxxx@mn1.bsc.es

Load the ruby module:

    $ module load ruby

Switch to the BOAST directory:

    $ cd boast

Compile BOAST:

    $ gem build BOAST.gemspec
    $ gem install --local --user-install BOAST-2.1.0.gem

Add ruby to your PATH (do it only once)

    $ echo 'export PATH=$PATH:/gpfs/home/bsc21/bsc21xxx/.gem/ruby/2.4.0/bin' >> ~/.bashrc

## Apply BOAST on nastin's split kernel

Here, we explain how to apply BOAST on one of Alya's kernels, `split`.
First, go into your Alya directory.

    $ cd /path/to/Alya

Go to the BOAST repository

    $ cd Boast/mod_nsi_assembly_split_oss/boast/Run

### Benchmark split on your platform

In order to determine the best code optimization and compilation options, you may run the following benchmark:

    $ ./benchmark.sh

It will generate pdf files that contain data about the best optimization and compilation options, both for the reference and the boast version. This can help you to select the boast version and the compilation options you want.

### Compile Alya with BOAST

First, you need to set up Alya's configuration file.

    $ cp /path/to/Alya/Executable/unix/configure.in/config_ifort.in /path/to/Alya/Executable/unix/config.in

It works for both ifortran and gfortran versions.

You need now to chose your optimization/compilation options among the following:

    Optimization:         -O2 or -O3 
    Unrolling:            -u if enabled, no flag if disabled
    Inlining:             -i [inlined/call/included]
    Vector size:          -v [number]

Then, launch the BOAST split kernel generation, adapting the options you want. Example:

    $ ./generate.sh -O3 -u -i inlined -v 4

This script will:

  - generate the corresponding fortran split source code,
  - copy it automatically in your Alya source directory,
  - modify your config.in file according to the optimization option and vector size you have defined. It will erase the exisiting optimization and vector size settings defined in this file.

Now, go to your Alya's `Executable/unix` directory and configure and compile Alya as usual.

    $ cd /path/to/Alya/Executable/unix
    $ ./configure -x parall nastin
    $ make metis4
    $ make -j4
