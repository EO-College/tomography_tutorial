[![Documentation Status][1]][2]

# EO-College tomography tutorial

This tutorial, developed by the [EO-College][3] learning initiative,
 explores Synthetic Aperture Radar (SAR) tomography 
with data from DLR's [F-SAR][4] system.
It consists of a Python package, containing several functions for processing and displaying the data, 
a Jupyter notebook and a test data set, which can be downloaded [here][5].  
Please follow the steps below to get started.

## Installation

The following subsections describe the installation process for different operating systems.
Please mind that this tutorial depends on Python 3.

#### Ubuntu

First we want to install [GDAL][9] to read our data. For this we add the [UbuntuGIS][10] package 
repository so we can install a more recent version than that supplied by Ubuntu.
After this we install GDAL together with its Python bindings:
```sh
sudo add-apt-repository ppa:ubuntugis/ppa
sudo apt-get update
sudo apt-get install gdal-bin python3-gdal
```

Next we install Tkinter for graphical support:
```sh
sudo apt-get install python3-tk
```

As a last step we install the tomography module including its direct Python package 
dependencies:

```sh
sudo python3 -m pip install tomography_tutorial
```

#### Windows

The easiest way to install Python and Jupyter on Windows is via [Anaconda][6]. 
Please make sure to install the Python 3 version.  
Once you have installed it, please add its installation directory to the PATH environment variable. 
See e.g. [here][7] for instructions.
Now we can install GDAL via Anaconda's own command line installation program:
```sh
conda install -c conda gdal
```

Finally we can install the tutorial package:
```sh
python -m pip install tomography_tutorial
```

## download of tutorial test data
Prior to starting the tutorial you need to download and unpack the data found 
[here][3].

## Starting the notebook

Now that everything is installed you can start the notebook via the tutorial Python module.
In the command prompt, start Python and execute the function [start][8]:
```Python
from tomography_tutorial import start
start('/your/custom/notebook.ipynb')
```
This will create a custom copy of the notebook if it does not yet exist and start it in the browser.
If the directory in which the custom notebook is to be stored does not yet exist, it is created 
automatically. Please mind that under Windows paths need to be separated with `\\` or `/`, 
a single backslash will cause an error.  
You now have a custom version of the tutorial, 
which you can modify as you like and restart later via function `start`.  
If you want to restore the original notebook, which was delivered with the Python package, just delete 
your custom version and run function `start` again.  

## API documentation

The documentation of the package functionality which is used in the notebook can be found 
[here][2].

[1]: https://readthedocs.org/projects/eocollege-tomography/badge/?version=latest
[2]: http://eocollege-tomography.readthedocs.io/en/latest/?badge=latest
[3]: https://eo-college.org/landingpage/
[4]: https://www.dlr.de/hr/en/desktopdefault.aspx/tabid-2326/3776_read-5691
[5]: https://eo-college.org/Data/Tomography/tomography_data.zip
[6]: https://conda.io/docs/user-guide/install/windows.html
[7]: https://www.computerhope.com/issues/ch000549.htm
[8]: https://eocollege-tomography.readthedocs.io/en/latest/tomography.html#tomography_tutorial.functions.start
[9]: https://www.gdal.org/
[10]: https://wiki.ubuntu.com/UbuntuGIS