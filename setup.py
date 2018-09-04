import os
import sys
from setuptools import setup, find_packages

directory = os.path.abspath(os.path.dirname(__file__))
if sys.version_info >= (3, 0):
    with open(os.path.join(directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
else:
    with open(os.path.join(directory, 'README.md')) as f:
        long_description = f.read()

setup(name='tomography_tutorial',
      packages=find_packages(),
      include_package_data=True,
      version='0.9',
      description='A tutorial package for Synthetic Aperture Radar Tomography',
      classifiers=[
          'Programming Language :: Python',
      ],
      python_requires='>3.0.0',
      install_requires=['numpy',
                        'jupyter',
                        'matplotlib',
                        'scipy'],
      url='https://github.com/EO-College/tomography_tutorial',
      author='John Truckenbrodt',
      author_email='john.truckenbrodt@uni-jena.de',
      license='MIT',
      zip_safe=False,
      long_description=long_description,
      long_description_content_type='text/markdown')
