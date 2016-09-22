.. _faq:

Frequently Asked Questions
==========================

When creating a photon list, sometimes pyXSIM just hangs and then after a while it crashes. What's happening?
-------------------------------------------------------------------------------------------------------------

This is likely to happen if you are simulating a thermal source and you have a very large effective area
and/or long exposure time. You may have to split the observation into a number of exposures, and then join them
together using :func:`~pyxsim.utils.merge_files`, or run the photon list generation in parallel on a larger machine.

What happens if my data source straddles a periodic boundary?
-------------------------------------------------------------

As of version 1.01, this is fine so long as your data source is a sphere or rectangular solid which is aligned with
the grid axes. The coordinates will be wrapped properly. For data sources with other shapes, this is currently not
implemented and a warning will be given to alert you of this issue.

