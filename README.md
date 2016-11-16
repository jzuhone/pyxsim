[![yt](http://img.shields.io/badge/powered%20by-yt-blue.svg?style=flat)](http://yt-project.org/)    [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)  
[![Build Status](https://travis-ci.org/jzuhone/pyxsim.svg?branch=master)](https://travis-ci.org/jzuhone/pyxsim)  
[![Coverage Status](https://coveralls.io/repos/github/jzuhone/pyxsim/badge.svg?branch=master)](https://coveralls.io/github/jzuhone/pyxsim?branch=master)  
# What is pyXSIM?

pyXSIM is a Python package for simulating X-ray observations from astrophysical sources.

X-rays probe the high-energy universe, from hot galaxy clusters to compact objects such as
neutron stars and black holes and many interesting sources in between. pyXSIM makes it
possible to generate synthetic X-ray observations of these sources from a wide variety of 
models, whether from grid-based simulation codes such as FLASH, Enzo, and Athena, to
particle-based codes such as Gadget and AREPO, and even from datasets that have been created
"by hand", such as from NumPy arrays. pyXSIM also provides facilities for manipulating the 
synthetic observations it produces in various ways, as well as ways to export the simulated
X-ray events to other software packages to simulate the end products of specific X-ray
observatories. 

# The Heritage of pyXSIM

pyXSIM is an implementation of the [PHOX](http://www.mpa-garching.mpg.de/~kdolag/Phox/)
algorithm, developed for constructing mock X-ray observations from SPH datasets by
Veronica Biffi and Klaus Dolag. There are two relevant papers:

[Biffi, V., Dolag, K., Bohringer, H., & Lemson, G. 2012, MNRAS, 420, 3545](http://adsabs.harvard.edu/abs/2012MNRAS.420.3545B)

[Biffi, V., Dolag, K., Bohringer, H. 2013, MNRAS, 428, 1395](http://adsabs.harvard.edu/abs/2013MNRAS.428.1395B)

pyXSIM had a previous life as the `photon_simulator` analysis module as a part of the
[yt Project](http://yt-project.org). pyXSIM still depends critically on yt to provide the
link between the simulation data and the algorithm for generating the X-ray photons. For
detailed information about the design of the algorithm in yt, check out
[the SciPy 2014 Proceedings](http://conference.scipy.org/proceedings/scipy2014/zuhone.html).

# Installing pyXSIM

pyXSIM can be installed in a few different ways. The simplest way is via the conda package if
you have the [Anaconda Python Distribution](https://store.continuum.io/cshop/anaconda/):

```
[~]$ conda install -c jzuhone pyxsim
```

The second way to install pyXSIM is via pip. pip will attempt to download the dependencies and
install them, if they are not already installed in your Python distribution:

```
[~]$ pip install pyxsim
```

Alternatively, to install into your Python distribution from [source](http://github.com/jzuhone/pyxsim):

```
[~]$ python setup.py install
```

# Getting Help

There are a number of ways to get help with pyXSIM. 

## Documentation

Documentation for pyXSIM lives at http://hea-www.cfa.harvard.edu/~jzuhone/pyxsim.

## Mailing List

There's a [pyXSIM Google Group](https://groups.google.com/forum/#!forum/pyxsim) to get help and 
discuss related matters.

## GitHub Issues Page

If you have a specific code issue that seems like a bug or have a feature or enhancement request,
the best place to note it is on the [GitHub issues page](http://github.com/jzuhone/pyxsim/issues)
so that we can keep track of it. 

## Slack Channel at the yt Project

The yt Project has a [Slack](http://www.slack.com) team which has a pyXSIM channel. Sign up for 
the [yt Slack team](http://yt-project.slack.com) to get onto the channel and ask questions. 
