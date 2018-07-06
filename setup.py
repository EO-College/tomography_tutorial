from setuptools import setup, find_packages
import os
import shutil

# Create .tomography in HOME - Directory
directory = os.path.join(os.path.expanduser('~'), '.tomography')

if not os.path.exists(directory):
    os.makedirs(directory)


setup(name='tomography',
      packages=find_packages(),
      include_package_data=True,
      version='0.1',
      description='A tutorial package for Synthetic Aperture Radar Tomography',
      classifiers=[
          'Programming Language :: Python',
      ],
      python_requires='>3.0.0',
      install_requires=['numpy',
                        'jupyter',
                        'matplotlib',
                        'scipy'],
      url='https://github.com/SAR-EDU/tomography.git',
      author='John Truckenbrodt',
      author_email='john.truckenbrodt@uni-jena.de',
      license='MIT',
      zip_safe=False)

# copy the tutorial notebook from the directory of the installed package to .tomography
shutil.copyfile('./tomography/tutorial.ipynb', os.path.join(directory, 'tutorial.ipynb'))
