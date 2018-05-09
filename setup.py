from setuptools import setup, find_packages

setup(name='BOFdat',
      version='0.1.3',
      download_url = 'https://github.com/peterldowns/mypackage/archive/0.3.tar.gz',
      description='Package to generate biomass objective function stoichiometric coefficients from experimental data',
      url='https://github.com/jclachance/BOFdat/',
      author='Jean-Christophe Lachance',
      author_email='jelachance@eng.ucsd.edu',
      license='MIT',
      packages=find_packages(),
      install_requires=['cobra>=0.8.0', 'BioPython','seaborn'])
