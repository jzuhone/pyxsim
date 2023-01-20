from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("pyxsim")
except PackageNotFoundError:
    # package is not installed
    pass


from pyxsim.event_list import EventList
from pyxsim.photon_list import (
    PhotonList,
    make_photons,
    project_photons,
    project_photons_allsky,
)
from pyxsim.source_models import (
    CIESourceModel,
    IGMSourceModel,
    LineSourceModel,
    NEISourceModel,
    PowerLawSourceModel,
)
from pyxsim.utils import (
    compute_elem_mass_fraction,
    compute_zsolar,
    create_metal_fields,
    merge_files,
)
