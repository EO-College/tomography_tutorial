from setuptools import setup, find_packages

setup(name='tomography',
      packages=find_packages(),
      include_package_data=True,
      version='0.1',
      description='a tutorial package for Synthetic Aperture Radar Tomography',
      classifiers=[
          'Programming Language :: Python',
      ],
      install_requires=['numpy',
                        'scipy'],
      url='https://github.com/johntruckenbrodt/tomography.git',
      author='John Truckenbrodt',
      author_email='john.truckenbrodt@uni-jena.de',
      license='MIT',
      zip_safe=False)
