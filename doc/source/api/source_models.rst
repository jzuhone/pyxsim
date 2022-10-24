Source Models API
=================

.. autoclass:: pyxsim.source_models.sources.SourceModel
    :members:
    :undoc-members:
    :exclude-members: setup_model, process_data, make_fluxf, set_pv, cleanup_model, compute_radius, setup_pbar

.. autoclass:: pyxsim.source_models.thermal_sources.CIESourceModel
    :members: make_spectrum

.. autoclass:: pyxsim.source_models.thermal_sources.NEISourceModel
    :members: make_spectrum

.. autoclass:: pyxsim.source_models.thermal_sources.IGMSourceModel
    :members: make_spectrum

.. autoclass:: pyxsim.source_models.power_law_sources.PowerLawSourceModel
    :members: make_spectrum
    :undoc-members:
    :exclude-members: setup_model, process_data, make_fluxf, cleanup_model

.. autoclass:: pyxsim.source_models.line_sources.LineSourceModel
    :members: make_spectrum
    :undoc-members:
    :exclude-members: setup_model, process_data, make_fluxf, cleanup_model

