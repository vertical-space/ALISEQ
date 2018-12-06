Instalation and setup
================================================================================



Download the source code
--------------------------------------------------------------------------------

.. code-block:: bash

   git clone https://github.com/vertical-space/ALISEQ
   git clone https://github.com/secastel/allelecounter.git ALISEQ/allelecounter

   export SRC_DIR=`pwd`
   conda env create -f ALISEQ/envs/aliseq.yaml -p ${SRC_DIR}/ALISEQ
   source activate aliseq

   export TEST_DIR=testdirectory # or whatever you want to call your working directory

   mkdir ${TEST_DIR}   
   cp ALISEQ/data/aliseq_params.txt ${TEST_DIR}
   cp ALISEQ/data/simulation_params.txt ${TEST_DIR}
   cp ALISEQ/data/targets.txt ${TEST_DIR}
   cd ${TEST_DIR}


Test if it's working
--------------------------------------------------------------------------------

.. code-block:: bash

   bash ../ALISEQ/aliseq.sh

Note: The first time you run this it will take a while to download the dependencies

