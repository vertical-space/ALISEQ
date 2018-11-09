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

   mkdir testdirectory   # or whatever you want to call your working directory

   cp ALISEQ/data/targets.txt testdirectory

   cp ALISEQ/data/targets.txt testdirectory

   cd testdirectory


Test if it's working
--------------------------------------------------------------------------------

.. code-block:: bash

   bash ../ALISEQ/aliseq.sh

Note: The first time you run this it will take a while to download the dependencies

