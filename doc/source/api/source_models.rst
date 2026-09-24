Source Models API
=================

.. autoclass:: pyxsim.source_models.sources.SourceModel
    :members:
    :undoc-members:
    :exclude-members: setup_model, process_data, make_fluxf, set_pv, compute_radius, setup_pbar

.. autoclass:: pyxsim.source_models.thermal_sources.collisional.CIESourceModel
    :members: make_spectrum

.. autoclass:: pyxsim.source_models.thermal_sources.collisional.NEISourceModel
    :members: make_spectrum

.. autoclass:: pyxsim.source_models.thermal_sources.photoionization.PionSourceModel
    :members: make_spectrum

.. autoclass:: pyxsim.source_models.power_law_sources.PowerLawSourceModel
    :members: make_spectrum
    :undoc-members:
    :exclude-members: setup_model, process_data, make_fluxf

.. autoclass:: pyxsim.source_models.line_sources.LineSourceModel
    :members: make_spectrum
    :undoc-members:
    :exclude-members: setup_model, process_data, make_fluxf
