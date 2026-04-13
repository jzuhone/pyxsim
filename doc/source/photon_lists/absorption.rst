.. _absorption:

Absorption Internal to a Source
===============================

.. note::

    This is a simple implementation of the effects of internal
    foreground absorption, with some caveats and approximations
    that are noted below. We hope to improve the fidelity of this
    feature in future releases.

It is trivial to include the effects of foreground Galactic
absorption from the Milky Way in mock observations of a 3D source
using pyXSIM, since we assume that there is one value of the
neutral hydrogen column density (in cm\ :sup:`-2`) in between
the observer and the source. Including the effects of absorption
internal to such a source, however, is more challenging. This
requires knowing the total column density of neutral hydrogen
along the line of sight between the observer and each point in
the source.

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

    normal = "x" # in this case, project along x-axis
    _, c = ds.find_min(("PartType1", "Potential")) # find potential minimum for center
    width = (1.0, "Mpc") # width of map
    depth = (1.0, "Mpc") # depth of map
    nwidth = 256 # number of cells along width
    ndepth = 256 # number of cells along depth
    outfile = "nH.h5" # file to write the map to
    field = ("gas", "H_p0_number_density") # the neutral hydrogen field, this is the default value

    pyxsim.make_column_density_map(
        ds,
        normal,
        c,
        width,
        depth,
        nwidth,
        ndepth,
        outfile,
        field=field,
    )

This will run through the selected region and create a cube of
maps of neutral hydrogen column density, one map for each point
in ``ndepth`` cells along the line of sight direction. The deeper
one goes into this cube, the farther away the source is, and the
column density will increase.

Once this file is created, it can be used in conjunction with
the :func:`~pyxsim.photon_list.project_photons` function to apply
absorption to the photons in the same manner as the foreground
Galactic absorption (and assuming the same models), but now including
the absorption from whatever neutral gas the photons encounter along
the sight line. This is handled by simply specifying the path to the
column density file that was just created to the optional keyword
argument ``column_file``:

.. code-block:: python

    pyxsim.project_photons(
        "my_photons",
        "my_events",
        "x",
        [30.0, 45.0],
        absorb_model="tbabs",
        nH=nH_sim,
        prng=prng,
        column_file="nH.h5",
    )

A worked example of how to include internal absorption in a
simulation of the CGM is shown in `More_Advanced_Thermal_Emission`.

Limitations of This Method
--------------------------

The end-user should note the following limitations of this
method, which we hope to rectify in a future release:

* This method does not take into account the velocity of the
  absorbing gas, which can Doppler-shift the absorption features.
* This method only uses the neutral hydrogen column density to
  determine the absorption under the assumption of an absorption
  model, and does not take into account the presence of other
  elements that can contribute to the absorption.
