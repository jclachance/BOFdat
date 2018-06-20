|Buildstatus| |License| |Documentation|

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

A full biomass objective function stoichiometric coefficients determination from experimental data fetched from literature for the *E.coli* model iML1515 is available in the Example folder. The files used are also provided. 

|Binder|

Documentation
~~~~~~~~~~~~~
The documentation and API for BOFdat is available on `Read the Docs`_ 

.. _Read the docs: http://BOFdat.readthedocs.org/


Cite
----

BOFdat: generating biomass objective function stoichiometric coefficients from experimental data

https://www.biorxiv.org/content/early/2018/01/05/243881


.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: https://github.com/jclachance/BOFdat/blob/master/LICENSE
.. |Documentation| image:: https://readthedocs.org/projects/BOFdat/badge/?version=master
    :target: https://bofdat.readthedocs.io/en/latest/index.html
.. |Buildstatus| image:: https://travis-ci.org/jclachance/BOFdat.svg?branch=master
    :target: https://travis-ci.org/jclachance/BOFdat

.. |Binder| image::https://mybinder.org/badge.svg
    :target: https://github.com/jclachance/BOFdat/blob/master/Example_usage/Example.ipynb


Author: Jean-Christophe Lachance
Date: 06-13-2018
Version: 0.1.4
