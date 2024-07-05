Installation
============

scTail is available through `PyPI <https://pypi.org/project/scTail/>`_.
To install, type the following command line and add ``-U`` for updates:

.. code-block:: bash

   pip install -U scTail

Alternatively, you can install from this GitHub repository for latest (often
development) version by the following command line:

.. code-block:: bash

   pip install -U git+https://github.com/StatBiomed/scTail

In either case, add --user if you don’t have the write permission for your Python environment.

If you run into any issues or errors are raised during the installation process, feel free to contact us at GitHub_.

.. _GitHub: https://github.com/StatBiomed/scTail/issues

In addition, you should also install paraclu_ and bedtools by youself to make sure scTail running smoothly.

.. _paraclu: https://gitlab.com/mcfrith/paraclu

.. code-block:: bash 

   conda install bioconda::paraclu
   conda install bioconda::bedtools

 
