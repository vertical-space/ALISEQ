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

   cd testdirectory

Test if it's working
--------------------------------------------------------------------------------

.. code-block:: bash

   echo 'ERR2353209' > targets.txt

   bash aliseq.sh

It can take quite a while to download all of the required third party software

