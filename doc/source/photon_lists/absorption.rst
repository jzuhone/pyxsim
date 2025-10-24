.. _absorption:

Absorption Internal to a Source
===============================

It is trivial to include the effects of foreground Galactic
absorption in mock observations of a 3D source using pyXSIM,
since we assume that there is one value of the neutral hydrogen
column density (in cm\ :sup:`-2`) for the entire source.

Including the effects of absorption internal to such a source,
however, is more challenging. This requires knowing the total
column density of neutral hydrogen along the line of sight between
the observer and each point in the source.

pyXSIM implements a simple method for including the effects of
internal absorption in a source for use when projecting a photon
list to an event list. The basic idea is that we construct a 3D
array of column densites, where the width and height correspond to
distances on the sky, and the depth corresponds to the distance
along the line of sight. The neutral hydrogen column density at
each point in the array is determined by projecting the density
along the line element from the observer to that point. This array
is then saved to disk. It can then be later used in conjunction
with the :meth:`~pyxsim.photon_list.project_photons` method, where
for each photon in the set there is a chance that it will be
absorbed based on the column density between its position in the
source and the observer.

The method for creating the column density array is the
:meth:`~pyxsim.internal_absorption.make_column_density_map` method,
which takes a dataset, normal vector, width, depth, and number
of cells in the width and depth directions to define the 3D array
of column densities:

.. code-block:: python
