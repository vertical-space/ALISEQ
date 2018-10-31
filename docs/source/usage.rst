Running the software
================================================================================

Edit the file targets.txt
--------------------------------------------------------------------------------

manually edit the file targets.txt, placing one accession number per line, e.g.

::

   ERR2352620
   ERR2352630
   ERR2352640
   ERR2352650
   ERR2352660
   ERR2352670

alternatively, if you have run the simulate reads pipeline, you can specify 
simulated data, e.g.

::

   simulation_1_x1
   simulation_1_x2

Making sure the names correspond to fastq files in the reads directory, e.g.

::

   simulation_1_x1.fastq.gz
   simulation_1_x2.fastq.gz
  
Then you can just relaunch the workflow with:

.. code-block:: bash

   bash aliseq.sh


