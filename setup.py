from setuptools import setup, find_packages

setup(name='BOFdat',
      version='0.1.6',
      download_url = 'https://github.com/peterldowns/mypackage/archive/0.1.6.tar.gz',
      description='Package to generate biomass objective function stoichiometric coefficients from experimental data',
      url='https://github.com/jclachance/BOFdat/',
      author='Jean-Christophe Lachance',
      author_email='jelachance@eng.ucsd.edu',
      license='MIT',
      packages=find_packages(),
      package_data={'BOFdat':['data/*.*']},
      include_package_data=True,
      install_requires=['cobra>=0.11.0','numpy>=1.13',
			'BioPython','seaborn',
			'scikit-learn>=0.18','scipy',
			'deap','pebble','matplotlib'])
