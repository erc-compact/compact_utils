# Mosaic for COMPACT: Multibeamformed observation simulation and interferometry characterization


This is a forked repo of the original code written by Weiwei Chen. Here is the original repo https://github.com/wchenastro/Mosaic/tree/master

The codes within this repo has been adapted for the COMPACT project.


Mosaic has the following dependencies
## Dependent

For python 3.8.5+

- numpy
- scipy
- matplotlib
- astropy
- nvector
- geographiclib
- katpoint

For python 2.7,  A docker instance is recommended, the content of Dockerfile list below:

```
FROM ubuntu:16.04

MAINTAINER Weiwei Chen wchen@mpifr-bonn.mpg.de

RUN apt-get update && \
    apt-get --no-install-recommends -y install \
    wget python-pip python-setuptools python-wheel \
    build-essential python-dev python-scipy python-numpy \
    python-matplotlib python-astropy

RUN pip install 'nvector==0.7.0' 'pillow==4.0.0' WCSAxes geographiclib katpoint
```



## Usage


```
python loop_grid_make_tiling.py
```

