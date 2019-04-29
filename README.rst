|License| |Documentation|

BOFdat
======
Generate biomass objective function for genome-scale models from experimental data.
BOFdat is a three step workflow that allows modellers to generate a complete biomass objective function *de novo* from experimental data:

1. Obtain stoichiometric coefficients for major macromolecules and calculate maintenance cost

2. Find coenzymes and inorganic ions

3. Find metabolic end goals


Significance
------------

Genome-scale metabolic models rely both on a defined media and a precise biomass objective function to generate reliable predictions of flux-states and gene essentiality. Generate a biomass objective that is specific to your organism of interest by incorporating experimental data and calculating stoichiometric coefficients. This package aims to produce an easy way to generate biomass stoichiometric coefficients that reflect experimental reality by incorporating weight fractions and relative abundances of macromolecules obtained from multiple OMICs datasets and finding specie-specific metabolic end goals. 

Installation
~~~~~~~~~~~~

Use pip to install BOFdat from `PyPi`_::

	pip install BOFdat


.. _PyPi: https://pypi.org/project/BOFdat/

Example use
~~~~~~~~~~~

A full biomass objective function stoichiometric coefficients determination from experimental data fetched from literature for the *E.coli* model *i*ML1515 is available in the Example folder. The files used are also provided. 


Documentation
~~~~~~~~~~~~~
The documentation and API for BOFdat is available on `Read the Docs`_ 

.. _Read the docs: http://BOFdat.readthedocs.org/


Cite
----

|BOFdat Generating biomass objective functions for genome-scale metabolic models from experimental data|_


.. _BOFdat Generating biomass objective functions for genome-scale metabolic models from experimental data: https://doi.org/10.1371/journal.pcbi.1006971
.. |BOFdat Generating biomass objective functions for genome-scale metabolic models from experimental data| replace:: BOFdat: Generating biomass objective functions for genome-scale metabolic models from experimental data

.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: https://github.com/jclachance/BOFdat/blob/master/LICENSE
.. |Documentation| image:: https://readthedocs.org/projects/BOFdat/badge/?version=master
    :target: https://bofdat.readthedocs.io/en/latest/index.html

Author: Jean-Christophe Lachance
Date: 06-13-2018
Version: 0.1.4
