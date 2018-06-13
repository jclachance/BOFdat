|Build_Status| |License| |Documentation|

BOFdat
======
Package to generate biomass objective function for genome-scale models from experimental data.


Significance
------------

Genome-scale metabolic models rely both on a defined media and a precise biomass objective function to generate reliable predictions of flux-states and gene essentiality. Generate a biomass objective that is specific to your organism of interest by incorporating experimental data and calculating stoichiometric coefficients. This package aims to produce an easy way to generate biomass stoichiometric coefficients that reflect experimental reality by incorporating weight fractions and relative abundances of macromolecules obtained from multiple OMICs datasets. 

Example use
-----------

A full biomass objective function stoichiometric coefficients determination from experimental data fetched from literature for the *E.coli* model *i*ML1515 is available in the Example folder. The files used are also provided. 

|Binder|]()]()

Documentation
-------------
Visit the [documentation](http://BOFdat.readthedocs.org/)



Cite
----
BOFdat: generating biomass objective function stoichiometric coefficients from experimental data

https://www.biorxiv.org/content/early/2018/01/05/243881


.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: https://github.com/jclachance/BOFdat/blob/master/LICENSE
.. |Documentation| image:: https://readthedocs.org/projects/BOFdat/badge/?version=master
    :target: https://bofdat.readthedocs.io/en/latest/index.html
.. |Binder| image::https://mybinder.org/badge.svg
    :target:https://mybinder.org/v2/gh/jclachance/BOFdat/blob/master/Example_usage/Example.ipynb


Author: Jean-Christophe Lachance
Date: 06-13-2018
Version: 0.1.3
