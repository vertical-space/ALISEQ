Running the read simulation pipeline
================================================================================


Modifying simulation parameters
--------------------------------------------------------------------------------

If desired, manually edit the file ``simulation_params.yaml``

Uncomment lines as required then adjust settings as desired. Use the default 
values (provided) as a starting point

It is necesary to set a new "output_tag" for each run, unless "clobber" is also 
used, in which case it will overwrite previous runfiles of the same name, if they 
exist.


Running the workflow
--------------------------------------------------------------------------------

launch the workflow with:

.. code-block:: bash

   bash ../ALISEQ/simulate_reads.sh


.. note::

   If you encounter the following error mesage (e.g. after a power failure or loss of network connection) ``Error: Directory cannot be locked...``, run the following cmd before retrying
   
   .. code-block:: bash

      snakemake --unlock -s ../ALISEQ/scripts/simulate_reads.snk

