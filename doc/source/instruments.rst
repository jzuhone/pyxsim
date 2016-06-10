.. _instruments:

Convolving Events with Instrumental Responses
=============================================

.. note::
-
    If you want to export simulated events to an instrument simulator, DO NOT
    convolve them with either an ARF or RMF from within ``pyxsim``. This will
    be the job of the external programs to do. 

MARX
----

As of MARX version 5.3.1, simulated events can be read in via ``SIMPUT`` files. 


SIMX
----

Events can be inputted into ``SIMX`` by use of a ``SIMPUT`` file. 