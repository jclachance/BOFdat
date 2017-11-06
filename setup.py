from setuptools import setup, find_packages

setup(name='BOFdat',
      version='0.3',
      description='Package to generate biomass objective function stoichiometric coefficients from experimental data',
      url='https://github.com/jclachance/BOFdat/',
      author='Jean-Christophe Lachance',
      author_email='jelachance@eng.ucsd.edu',
      license='MIT',
      packages=find_packages(),
      zip_safe=False, install_requires=['cobra'])
