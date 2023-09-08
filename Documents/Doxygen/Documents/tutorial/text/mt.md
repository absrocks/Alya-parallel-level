# Running Alya on Minotauro {#mt}

## Compilation with PGI

Before trying to compile Alya, you have to load the following modules:

    module load cuda/8.0 pgi/17.10

And set this environment variable to use the pgfortran compiler through the mpi wrapper.

    export METIS_CC=pgcc

Type these commands for each new session on Marenostrum, before trying to compile Alya.

The configurations files that are used in the following sections are available in the `Executables/unix/configure.in/mt/` directory.

### Compiling Alya for CPU

Let's create a CPU Executable directory. In the Alya directory:

    cd Executables
    cp -r unix unix-cpu
    cd unix-cpu

Now, let's configure Alya. Copy the configuration file. Edit it if necessary.
    cp configure.in/mt/config-cpu.in config.in

Run configure with the right modules:
    ./configure -x parall nastin

Compile metis4

    make metis4

Compile Alya

    make

### Compiling Alya for GPU

Let's create a GPU Executable directory. In the Alya directory:

    cd Executables
    cp -r unix unix-gpu
    cd unix-gpu

Now, let's configure Alya. Copy the configuration file. Edit it if necessary.
    cp configure.in/mt/config-gpu.in config.in

Run configure with the right modules:

    ./configure -x parall nastin

You don't need to recompile metis.

Compile Alya

    make -j16

## Configure your case

In the `.dat` file, you will have to specify the number of cores per GPUS. In the case of Minotauro, it is 4 for the clusters k80 and 6 for the clusters m2090.

    PARALL_SERVICE:         On
      PARTITION_TYPE:       FACES
      PARTITIONING:
        METHOD:             SFC
      END_PARTITIONING
      CPGPU:                4
    END_PARALL_SERVICE 
    
In the `nsi.dat` file, you'll have to set the convective term to EMAC:

    PROBLEM_DEFINITION
      [...]
      CONVECTIVE_TERM:          EMAC
      [...]
    END_PROBLEM_DEFINITION

And set the following values in the numerical treatment section:

    NUMERICAL_TREATMENT
      [...]
      DIRICHLET: ALGORITHM
      GRAD_DIV= ON
      ASSEMBLY:   GPU
      [...]
    END_NUMERICAL_TREATMENT

## Configure the job

As we have two versions of Alya, we must provide mpirun with the information about which process runs the GPU or the CPU version.
Here, an example of such a job script allocating two nodes.

### submit.sh

This is the main script that you submit with `sbatch submit.sh`

    #!/bin/bash
    #SBATCH --job-name=sphere
    #SBATCH --output=sphere_%j.out
    #SBATCH --error=sphere_%j.err
    #SBATCH --ntasks=32
    #SBATCH --time=01:00:00
    #SBATCH --cpus-per-task=1
    #SBATCH --ntasks-per-node=16
    #SBATCH --gres gpu:4
    #SBATCH --constraint=k80

    export ROMIO_HINTS=./io_hints
    export I_MPI_EXTRA_FILESYSTEM_LIST=gpfs
    export I_MPI_EXTRA_FILESYSTEM=on
    module load cuda/8.0 pgi/17.10
    export NO_STOP_MESSAGE=1

    srun --multi-prog input.conf

### input.conf

This file is a wrapper that enables to retrieve for each MPI process the task number using the variable `%t`.

    * ./Alya.sh %t

### Alya.sh

This script determines which version of Alya each process will launch, regarding its task number.
It only works well if the nodes are exclusive, and with the assumption that the task numbers are sorted by sockets.

    #!/bin/bash
    #Parameters
    case=sphere                             #Case name
    only_cpu="false"                        #Only run with CPU
    only_gpu="false"                        #Only run with GPU
    threads_per_gpu=1                       #MPI Thread by GPU
    cpu_per_gpu=4                           #CPU count per GPU
    alya_dir=~/Alya/Executables/              #Alya Executable Directory
    #-----------Beside this section, should not be modified-----------
    cpu_dir=$alya_dir/unix-cpu
    gpu_dir=$alya_dir/unix-gpu
    exec=Alya.x
    cpu=$cpu_dir/$exec
    gpu=$gpu_dir/$exec
    ngpu=$threads_per_gpu
    ncpu=$(($cpu_per_gpu-$ngpu))
    TASK=$1
    if [ "$only_cpu" == "true"  ] || [ "$only_gpu" == "false" ] && [ `echo "$TASK % $cpu_per_gpu" | bc` -lt $ncpu ]
    then
      echo $TASK is running on CPU
      $cpu $case
    else
      echo $TASK is running on GPU
      $gpu $case
    fi


