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







Other stuff
--------------------------------------------------------------------------------

this is a demo of how to generate a document in sphinx

.. function:: secondFunc(arg1)

    not a real function!

.. this is an rst comment. it explains what the code does, but does not 
   show up in the html document.

   it can continue with indentation.

   and so on.

.. The following syntax imports a module and documents all of its members by using their docstrings

.. code-block:: python

    import random

    random.seed(42)

    def randomGene(m=300,M=10000):
        """create a random ATCG string with length randomly chosen from within the range [m,M]
        """
        length = 0
        while not m < length < (M+1):
	    length = int(np.random.normal(3392,2600)) # mu, sigm
        return ''.join([random.choice('ATCG') for i in range(length)])

    randomGene()

.. note::

    this is an example of a note. i'm not really sure what it does yet, so this is also a test of that.

    this is paragraph 2 of my fancy note.

.. warning::

    I'm pretty sure i know what this does!

.. class:: Request

    this is one way to document a class

    .. method:: func1()

    description of the func1 method

    .. attribute:: attr1

    description of the attr1 attribute

.. _Link: https://vertical-space.github.io/ALISEQ/


running the simulation
--------------------------------------------------------------------------------

.. figure:: http://fossilshelf.com/images/museum/IMG_7713.jpg
   :alt: this is an insect

   This is the caption of the figure (a simple paragraph).


running the ASE workflow
--------------------------------------------------------------------------------

* list item 1
* list item 2
* list item 3




