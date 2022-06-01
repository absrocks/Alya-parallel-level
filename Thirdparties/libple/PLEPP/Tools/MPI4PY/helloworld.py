#!/usr/bin/env python
"""
Parallel Hello World
"""

import sys
sys.path.append("/home/bsc21/bsc21704/z2015/REPOSITORY/MPI4PY131/Execs/lib64/python/") 

from mpi4py import MPI
import sys

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

sys.stdout.write(
    "Hello, World! I am process %d of %d on %s.\n"
    % (rank, size, name))
