Running the allele specific expression (ASE) pipeline
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

Make sure the names correspond to the prefix of paired fastq files in the reads directory, e.g.

::

   reads/simulation_1_x1.1.fastq.gz
   reads/simulation_1_x1.2.fastq.gz
   reads/simulation_1_x2.1.fastq.gz
   reads/simulation_1_x2.2.fastq.gz
  
Then launch the workflow with:

.. code-block:: bash

   bash ../ALISEQ/aliseq.sh


