[![Documentation Status](https://readthedocs.org/projects/eocollege-tomography/badge/?version=latest)](http://eocollege-tomography.readthedocs.io/en/latest/?badge=latest)
 
# tomography

A tutorial for Synthetic Aperture Radar Tomography  

This tutorial is still in the making as part of the 
<a href="https://eo-college.org/landingpage/" target="_blank" rel="noopener noreferrer">EO College</a> learning initiative.  
Please stay tuned..
## Installation

The following subsections descripbe the installation process for different operating systems.
Please mind that this tutorial depends on Python 3.

#### Ubuntu

First we want to install GDAL to read our data. For this we add the ubuntugis package 
repository so we can install a more recent version than that supplied by Ubuntu.  
After this we install GDAL together with its Python bindings.
```sh
sudo add-apt-repository ppa:ubuntugis/ppa
sudo apt-get update
sudo apt-get install gdal-bin python3-gdal
```

Next we install Tkinter for graphical support and git for package version control:
```sh
sudo apt-get install python3-tk
sudo apt-get install git
```

As a last step we install the tomography module including its direct Python package 
dependencies:

```sh
sudo python3 -m pip install git+https://github.com/SAR-EDU/tomography.git
```

#### Windows

The easiest way to install Python and Jupyter on Windows is via 
<a href="https://conda.io/docs/user-guide/install/windows.html" target="_blank" rel="noopener noreferrer">Anaconda</a>. 
Please make sure to install the Python 3 version.  
Once you have installed it, add its installation directory to 
the PATH environment variable. See e.g. 
<a href="https://www.computerhope.com/issues/ch000549.htm" target="_blank" rel="noopener noreferrer">here</a> for instructions.
We further need the versioning system git, which can be downloaded from 
<a href="https://git-scm.com/downloads" target="_blank" rel="noopener noreferrer">here</a>.  
Now we can install GDAL via Anaconda's own command line installation program:
```sh
conda install -c conda gdal
```

Finally we can install the tutorial package:
```sh
pip install git+https://github.com/SAR-EDU/tomography.git
```
