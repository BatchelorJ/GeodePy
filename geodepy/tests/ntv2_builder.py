"""
Helper for assembling small, byte-accurate synthetic NTv2 (.gsb) files for
tests, so geodepy.ntv2reader can be exercised without needing a real
(tens-to-hundreds-of-MB, separately licensed) grid on disk.

Mirrors the binary layout read by geodepy.ntv2reader.read_ntv2_file: a
176-byte overview header, one 176-byte header per subgrid, followed by that
subgrid's gs_count grid-shift node records (4 x little-endian float32 each).
"""

import os
import struct
from dataclasses import dataclass
from typing import Dict


GRIDS_DIR_ENV_VAR = "GEODEPY_NTV2_GRIDS_DIR"


@dataclass
class Ntv2ControlPoint:
    """
    An independently verified lat/lon -> NTv2 field value(s) check against
    a real grid file on disk. See
    geodepy/tests/ntv2_real_world_control_points.py (committed, organised
    into sections per grid product) and
    ntv2_real_world_control_points_local.py (optional, gitignored, for
    your own ad hoc grids) for where these are registered, and
    TestNtv2RealWorldControlPoints in test_ntv2reader.py for how they're
    run.

    :param grid_path: either a bare filename (e.g. "A66_National_13_09_01.gsb"),
        resolved at test time against the GEODEPY_NTV2_GRIDS_DIR
        environment variable - use this for grids anyone can download from
        a shared source like https://github.com/icsm-au/transformation_grids
        so the entry stays portable and committable - or an absolute path,
        used as-is, for a grid only you have (keep those in the local
        overlay file, not the committed one)
    :param lat: latitude of the control point (decimal degrees)
    :param lon: longitude of the control point (decimal degrees)
    :param expected: dict of any subset of "field_1".."field_4" to their
        independently verified values - i.e. NOT values obtained by
        running interpolate_ntv2() yourself, or the test can't catch a
        regression in interpolate_ntv2
    :param source: short citation of where `expected` came from
    :param label: optional short name; derived from grid_path if omitted
    :param method: "bicubic" (default) or "bilinear"
    :param places: decimal places of tolerance for assertAlmostEqual
    """

    grid_path: str
    lat: float
    lon: float
    expected: Dict[str, float]
    source: str
    label: str = ""
    method: str = "bicubic"
    places: int = 4

    def __post_init__(self):
        if not self.label:
            self.label = os.path.splitext(os.path.basename(self.grid_path))[0]


def resolve_grid_path(grid_path):
    """
    Resolve a Ntv2ControlPoint.grid_path to an actual file path.

    An absolute path, or one that already contains a directory component,
    is returned unchanged. A bare filename is looked up under the
    GEODEPY_NTV2_GRIDS_DIR environment variable instead, so committed
    control points never need to name a machine-specific location.

    :return: the resolved path, or None if grid_path is a bare filename
        and GEODEPY_NTV2_GRIDS_DIR isn't set.
    """
    if os.path.isabs(grid_path) or os.path.dirname(grid_path):
        return grid_path
    grids_dir = os.environ.get(GRIDS_DIR_ENV_VAR)
    if not grids_dir:
        return None
    return os.path.join(grids_dir, grid_path)


def _str_record(name, value, value_len=8):
    name_bytes = name.encode("ascii").ljust(8, b"\x00")[:8]
    value_bytes = value.encode("ascii").ljust(value_len, b"\x00")[:value_len]
    return name_bytes + value_bytes


def _int_record(name, value):
    name_bytes = name.encode("ascii").ljust(8, b"\x00")[:8]
    value_bytes = value.to_bytes(4, byteorder="little") + b"\x00" * 4
    return name_bytes + value_bytes


def _double_record(name, value):
    name_bytes = name.encode("ascii").ljust(8, b"\x00")[:8]
    value_bytes = struct.pack("d", value)
    return name_bytes + value_bytes


def build_overview_header(
    num_file,
    num_orec=11,
    num_srec=11,
    gs_type="SECONDS",
    version="NTv2.0",
    system_f="GDA94",
    system_t="GDA2020",
    major_f=6378137.0,
    minor_f=6356752.314140,
    major_t=6378137.0,
    minor_t=6356752.314140,
):
    return b"".join(
        [
            _int_record("NUM_OREC", num_orec),
            _int_record("NUM_SREC", num_srec),
            _int_record("NUM_FILE", num_file),
            _str_record("GS_TYPE", gs_type),
            _str_record("VERSION", version),
            _str_record("SYSTEM_F", system_f),
            _str_record("SYSTEM_T", system_t),
            _double_record("MAJOR_F", major_f),
            _double_record("MINOR_F", minor_f),
            _double_record("MAJOR_T", major_t),
            _double_record("MINOR_T", minor_t),
        ]
    )


def build_subgrid_header(
    sub_name,
    s_lat,
    n_lat,
    e_long,
    w_long,
    lat_inc,
    long_inc,
    gs_count,
    parent="NONE",
    created="01011990",
    updated="01011990",
):
    return b"".join(
        [
            _str_record("SUB_NAME", sub_name),
            _str_record("PARENT", parent),
            _str_record("CREATED", created),
            _str_record("UPDATED", updated),
            _double_record("S_LAT", s_lat),
            _double_record("N_LAT", n_lat),
            _double_record("E_LONG", e_long),
            _double_record("W_LONG", w_long),
            _double_record("LAT_INC", lat_inc),
            _double_record("LONG_INC", long_inc),
            _int_record("GS_COUNT", gs_count),
        ]
    )


def build_node_records(nodes):
    """
    :param nodes: rows ordered south-to-north (s_lat -> n_lat), each a list
        of (field_1, field_2, field_3, field_4) tuples ordered
        east-to-west as stored on disk (i.e. e_long -> w_long).
    """
    out = bytearray()
    for row in nodes:
        for field_1, field_2, field_3, field_4 in row:
            out += struct.pack("f", field_1)
            out += struct.pack("f", field_2)
            out += struct.pack("f", field_3)
            out += struct.pack("f", field_4)
    return bytes(out)


class Ntv2SubgridSpec:
    def __init__(
        self,
        sub_name,
        s_lat,
        n_lat,
        e_long,
        w_long,
        lat_inc,
        long_inc,
        nodes,
        parent="NONE",
        created="01011990",
        updated="01011990",
    ):
        self.sub_name = sub_name
        self.s_lat = s_lat
        self.n_lat = n_lat
        self.e_long = e_long
        self.w_long = w_long
        self.lat_inc = lat_inc
        self.long_inc = long_inc
        self.nodes = nodes
        self.parent = parent
        self.created = created
        self.updated = updated

    @property
    def gs_count(self):
        return sum(len(row) for row in self.nodes)


def build_ntv2_bytes(subgrids, **overview_kwargs):
    """
    :param subgrids: list of Ntv2SubgridSpec, in the order they should
        appear in the file
    :param overview_kwargs: forwarded to build_overview_header (e.g.
        gs_type, version, system_f, system_t, major/minor axes)
    :return: bytes of a complete, minimal NTv2 (.gsb) file
    """
    body = bytearray(build_overview_header(num_file=len(subgrids), **overview_kwargs))
    for sg in subgrids:
        body += build_subgrid_header(
            sub_name=sg.sub_name,
            s_lat=sg.s_lat,
            n_lat=sg.n_lat,
            e_long=sg.e_long,
            w_long=sg.w_long,
            lat_inc=sg.lat_inc,
            long_inc=sg.long_inc,
            gs_count=sg.gs_count,
            parent=sg.parent,
            created=sg.created,
            updated=sg.updated,
        )
        body += build_node_records(sg.nodes)
    return bytes(body)


def write_ntv2_file(path, subgrids, **overview_kwargs):
    data = build_ntv2_bytes(subgrids, **overview_kwargs)
    with open(path, "wb") as f:
        f.write(data)
    return path
