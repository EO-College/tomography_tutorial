[![Documentation Status](https://readthedocs.org/projects/eocollege-tomography/badge/?version=latest)](http://eocollege-tomography.readthedocs.io/en/latest/?badge=latest)
 
# tomography-tutorial

A tutorial for Synthetic Aperture Radar Tomography  

This tutorial is still in the making as part of the 
[EO College](https://eo-college.org/landingpage/) learning initiative.  
Please stay tuned..
## Installation

The following subsections describe the installation process for different operating systems.
Please mind that this tutorial depends on Python 3.

#### Ubuntu

First we want to install GDAL to read our data. For this we add the ubuntugis package 
repository so we can install a more recent version than that supplied by Ubuntu.
After this we install GDAL together with its Python bindings:
```sh
sudo add-apt-repository ppa:ubuntugis/ppa
sudo apt-get update
sudo apt-get install gdal-bin python3-gdal
```

Next we install Tkinter for graphical support and git for package version control:
```sh
sudo apt-get install python3-tk git
```

As a last step we install the tomography module including its direct Python package 
dependencies:

```sh
sudo python3 -m pip install git+https://github.com/SAR-EDU/tomography.git
```

#### Windows

The easiest way to install Python and Jupyter on Windows is via 
[Anaconda](https://conda.io/docs/user-guide/install/windows.html). 
Please make sure to install the Python 3 version.  
Once you have installed it, please add its installation directory to the PATH environment variable. 
See e.g. [here](https://www.computerhope.com/issues/ch000549.htm) for instructions.
Now we can install GDAL and git via Anaconda's own command line installation program:
```sh
conda install -c conda gdal git
```

Finally we can install the tutorial package:
```sh
pip install git+https://github.com/SAR-EDU/tomography.git
```
## Starting the notebook

Now that everything is installed you can start the notebook via the tutorial Python module.
In the command prompt, start Python and execute the function `start`:
```Python
from tomography_tutorial import start
start('/your/custom/notebook.ipynb')
```
This will create a custom copy of the notebook if it does not exist and start it in the browser.
If the directory in which the custom notebook is to be stored does not yet exist, it is created 
automatically. Please mind that under Windows paths need to be separated with `\\` or `/`.  
A single backslash will cause an error.  
You now have a custom version of the tutorial, 
which you can modify as you like and restart via function `start`.  
If you want to restore the original notebook, which was delivered with the Python package, just delete 
your custom version and run function `start` again.
