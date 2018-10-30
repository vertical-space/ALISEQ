Instructions for use
================================================================================


Editing the config file
--------------------------------------------------------------------------------

this is a demo of how to generate a document in sphinx

.. function:: secondFunc(arg1)

    not a real function!

.. this is an rst comment. it explain swhat the code here does, but should not show up in the actual document.

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

.. _Link: https://github.com/vertical-space/practise-public



running the simulation
--------------------------------------------------------------------------------

.. image:: http://fossilshelf.com/images/museum/IMG_7713.jpg
   :alt: this is an insect

.. figure:: ../media/IMG_6774.JPG
   :alt: this is an insect

   This is the caption of the figure (a simple paragraph).




running the ASE workflow
--------------------------------------------------------------------------------

* list item 1
* list item 2
* list item 3




