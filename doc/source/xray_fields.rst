.. _xray-fields:

X-ray Fields for yt
===================

The simplest use of a :class:`~pyxsim.source_models.sources.SourceModel` in pyXSIM is 
to create fields of X-ray emission that can be used in yt to compute emissivities, 
luminosities, or intensities. These can be computed for geometric objects such as 
spheres or boxes or can be projected along a sight line. There are two methods to do
this which are described below.

Source Fields
-------------

One may want to compute fields for the emissivity or the luminosity of a particular
source in its rest frame. This can be done using the 
:meth:`~pyxsim.source_models.sources.SourceModel.make_source_fields` method. This
method takes a :class:`~yt.data_objects.static_output.Dataset` object from yt, 
and the minimum and maximum energies of the band you want to create fields for. 
The example below shows creating source fields for thermal emission from a simulation
of a galaxy cluster: 

.. code-block:: python

    import yt
    import pyxsim
    
    ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100", default_species_fields="ionized")
    
    emin = 0.1
    emax = 10.0
    nbins = 1000
    Zmet = 0.3 # this dataset does not have a metallicity field, so assume 0.3 Zsolar
    source_model = pyxsim.CIESourceModel("apec", emin, emax, nbins, Zmet)
    
    # arguments are the Dataset, and the emin and emax of the 
    # band 
    xray_fields = source_model.make_source_fields(ds, 0.5, 7.0)

The fields are created for the :class:`~yt.data_objects.static_output.Dataset`
``ds``, and their names are returned in the ``xray_fields`` list:

.. code-block:: python

    print(xray_fields)

.. code-block:: pycon

    [('gas', 'xray_emissivity_0.5_7.0_keV'), 
     ('gas', 'xray_luminosity_0.5_7.0_keV'), 
     ('gas', 'xray_photon_emissivity_0.5_7.0_keV')]
    
Three fields have been created--one for the X-ray emissivity in the chosen band in
:math:`\rm{erg}~\rm{cm}^{-3}~\rm{s}^{-1}`, another for the X-ray luminosity in the
chosen band in :math:`\rm{erg}~\rm{s}^{-1}`, and another for the X-ray photon
emissivity in :math:`\rm{photon}~\rm{cm}^{-3}~\rm{s}^{-1}`. These fields exist in
the same way as any other field in yt, and can be used in the same ways. 

Querying emissivity values in a sphere:

.. code-block:: python

    sp = ds.sphere("c", (500.0, "kpc"))
    print(sp['gas', 'xray_emissivity_0.5_7.0_keV'])

.. code-block:: pycon

    [6.75018212e-30 6.63582106e-30 6.45686636e-30 ... 2.59468150e-30
     2.55886161e-30 2.65063999e-30] erg/(cm**3*s)

Summing luminosity in a sphere:

.. code-block:: python

    print(sp.sum(("gas", "xray_luminosity_0.5_7.0_keV")))

.. code-block:: pycon

    unyt_quantity(7.73753352e+44, 'erg/s')

Projecting the photon emissivity along a sight line:

.. code-block:: python

    prj = yt.ProjectionPlot(ds, "z", xray_fields[-1], width=(0.5, "Mpc"))
    prj.save()

.. image:: _images/projected_emiss.png