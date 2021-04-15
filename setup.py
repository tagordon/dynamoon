from setuptools import setup 

import pathlib
HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(name='dynamoon', 
      version='0.1',
      description='Photodynamics for exoplanet+moon systems',
      long_description=README,
      long_description_content_type="text/markdown",
      url='http://github.com/tagordon/dynamoon',
      author='Tyler Gordon',
      author_email='tagordon@uw.edu', 
      license='MIT',
      packages=['dynamoon'],
      install_requires=['numpy',
                        'scipy',
                        'astropy'],
      zip_safe=False)