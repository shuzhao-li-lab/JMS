from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
  name='jms-metabolite-services',
  version='0.2.0',

  author='Shuzhao Li',
  author_email='shuzhao.li@gmail.com',
  description='conversion, search of metabolic models and metabolomics data',
  long_description_content_type="text/markdown",
  long_description=long_description,
  url='https://github.com/shuzhao-li/JMS',
  license='Python',

  keywords='metabolomics, chemistry, bioinformatics, mass spectrometry',

  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3.7',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
  ],

  # changed from earlier setuptools
  packages=find_packages(
    include=['*', '']
  ),
  include_package_data=True,
  install_requires=[
    'mass2chem',
  ],

  python_requires='>=3',

)