- The script is updated for Jureca supercomputer. 

- Jube help is in http://apps.fz-juelich.de/jsc/jube/jube2/docu/index.html

- typical jube performance evaluation run command : 
jube -vvv run AlyaPerfEval.xml --tag ref mem darshan scatter compact scalasca no-vec no-fma

- known issue : PAPI run has way too much overhead, despite a big filtering list. needs more filtering to be relevant.
