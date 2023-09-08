# Tracing Alya with Extrae {#extrae}

## Instrumentation

If you do not instrumente the code, the whole execution will be traced by default.
It can be a problem if you run your simulation during a large time, because the trace generated will be voluminous and difficult to analyze.

To reduce the trace size, you can focus on a specific code part.

To do this, you need:

1) to load the extrae module in the fortran file to instrument:

    #IFDEF EXTRAE
      use extrae_module
    ENDIF

2) to (re)start extrae when entering to the code part to trace

    #ifdef EXTRAE
      call extrae_restart()
    #endif

3) you may trace one or several user events:

    call extrae_eventandcounters(type,int(value,8)) 
    
  with type the event type and value its value (both integers).

4) to stop extrae once the code part to trace ends.
    #ifdef EXTRAE
    call extrae_shutdown()
    #endif

## Configuration

To compile Alya with Extrae, you need to edit your configuration file by adding this information:

    ###################################################################
    #                          EXTRAE FLAGS                           #
    ###################################################################
    # Uncomment the following line to compile Alya using extrae
    # Compiler used to compile extrae module (make extrae)
    EXTRAE_FC=ifort
    # Extrae installation directory (for linking) (not necessary if loading extrae using module load extrae)
    #EXTRAE_HOME=/path/to/extrae
    # Extrae sources directory (for compiling Extrae module)
    EXTRAE_SRC=/path/to/sources/extrae/
    #@Linking with Extrae
    #EXTRALIB   := $(EXTRALIB) -L${EXTRAE_HOME}/lib/ -lmpitracef extrae_module.o
    #@Enabling Extrae user calls (normal)
    #CSALYA   := $(CSALYA) -DEXTRAE

Then, execute configure

    ./configure -x nastin parall

## Compilation

Compile the extrae module:

    make extrae

Compile Alya

    make


