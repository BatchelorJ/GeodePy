import io
import os
import re
import struct
import tempfile
import unittest
from importlib import import_module

from geodepy.ntv2reader import (
    bicubic_interpolation,
    bilinear_interpolation,
    interpolate_ntv2,
    read_node,
    read_ntv2_file,
)
from geodepy.tests.ntv2_builder import (
    GRIDS_DIR_ENV_VAR,
    Ntv2SubgridSpec,
    resolve_grid_path,
    write_ntv2_file,
)


class TestBilinearInterpolation(unittest.TestCase):
    def test_corners_return_node_values(self):
        n1, n2, n3, n4 = 11.0, 22.0, 33.0, 44.0
        self.assertAlmostEqual(bilinear_interpolation(n1, n2, n3, n4, 0, 0), n1)
        self.assertAlmostEqual(bilinear_interpolation(n1, n2, n3, n4, 1, 0), n2)
        self.assertAlmostEqual(bilinear_interpolation(n1, n2, n3, n4, 0, 1), n3)
        self.assertAlmostEqual(bilinear_interpolation(n1, n2, n3, n4, 1, 1), n4)

    def test_midpoint_matches_hand_calculation(self):
        # a0=0, a1=10, a2=20, a3=10 -> p = 0 + 10*x + 20*y + 10*x*y
        result = bilinear_interpolation(0, 10, 20, 40, 0.5, 0.5)
        self.assertAlmostEqual(result, 17.5)


class TestBicubicInterpolation(unittest.TestCase):
    # These three corners are derived analytically from the module's fixed
    # inverse-coefficient matrix (cinv): at each of (x, y) below, every term
    # involving the derivative/cross-derivative inputs (n5..n16) cancels
    # out of the sum, so the patch must reproduce the raw node value
    # exactly no matter what the surrounding nodes are.
    def setUp(self):
        self.nodes = {i: float(i) * 1.1 for i in range(1, 17)}

    def _call(self, x, y):
        n = self.nodes
        return bicubic_interpolation(
            n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8],
            n[9], n[10], n[11], n[12], n[13], n[14], n[15], n[16],
            x, y,
        )

    def test_corner_x0_y0_returns_n1(self):
        self.assertAlmostEqual(self._call(0, 0), self.nodes[1])

    def test_corner_x1_y0_returns_n2(self):
        self.assertAlmostEqual(self._call(1, 0), self.nodes[2])

    def test_corner_x0_y1_returns_n4(self):
        self.assertAlmostEqual(self._call(0, 1), self.nodes[4])


class TestReadNode(unittest.TestCase):
    def test_reads_four_floats_in_order(self):
        values = (1.5, -2.25, 0.0, 100.0)
        data = b"".join(struct.pack("f", v) for v in values)
        result = read_node(io.BytesIO(data))
        for expected, actual in zip(values, result):
            self.assertAlmostEqual(expected, actual, places=5)


def _node_value(row, col):
    base = row * 10 + col
    return (float(base), float(-base), 1.0, 2.0)


class TestReadAndInterpolateSingleSubgrid(unittest.TestCase):
    """
    Exercises read_ntv2_file and interpolate_ntv2 against a minimal,
    synthetic 4x4-node grid (one degree spacing, covering 0-3 deg lat/lon)
    so they can be verified without a real-world grid file.
    """

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmpdir.name, "synthetic.gsb")
        self.nodes = [[_node_value(row, col) for col in range(4)] for row in range(4)]
        self.subgrid_spec = Ntv2SubgridSpec(
            sub_name="TESTGRID",
            s_lat=0.0,
            n_lat=3 * 3600.0,
            e_long=0.0,
            w_long=3 * 3600.0,
            lat_inc=3600.0,
            long_inc=3600.0,
            nodes=self.nodes,
            created="15031995",
            updated="20062001",
        )
        write_ntv2_file(
            self.path,
            [self.subgrid_spec],
            gs_type="SECONDS",
            version="NTv2.0",
            system_f="GDA94",
            system_t="GDA2020",
            major_f=6378137.0,
            minor_f=6356752.314140,
            major_t=6378137.0,
            minor_t=6356752.314140,
        )
        self.grid = read_ntv2_file(self.path)

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_overview_header_fields(self):
        g = self.grid
        self.assertEqual(g.num_orec, 11)
        self.assertEqual(g.num_srec, 11)
        self.assertEqual(g.num_file, 1)
        self.assertEqual(g.gs_type, "SECONDS")
        self.assertEqual(g.version, "NTv2.0")
        self.assertEqual(g.system_f, "GDA94")
        self.assertEqual(g.system_t, "GDA2020")
        self.assertAlmostEqual(g.major_f, 6378137.0)
        self.assertAlmostEqual(g.minor_f, 6356752.314140)
        self.assertAlmostEqual(g.major_t, 6378137.0)
        self.assertAlmostEqual(g.minor_t, 6356752.314140)

    def test_subgrid_header_fields(self):
        sg = self.grid.subgrids["TESTGRID"]
        self.assertEqual(sg.parent, "NONE")
        self.assertEqual(sg.created, "15/03/1995")
        self.assertEqual(sg.updated, "20/06/2001")
        self.assertAlmostEqual(sg.s_lat, 0.0)
        self.assertAlmostEqual(sg.n_lat, 3 * 3600.0)
        self.assertAlmostEqual(sg.e_long, 0.0)
        self.assertAlmostEqual(sg.w_long, 3 * 3600.0)
        self.assertAlmostEqual(sg.lat_inc, 3600.0)
        self.assertAlmostEqual(sg.long_inc, 3600.0)
        self.assertEqual(sg.gs_count, 16)

    def test_interpolate_on_node_bilinear(self):
        # (1 deg, -1 deg) lands exactly on node (row=1, col=1)
        result = interpolate_ntv2(self.grid, 1.0, -1.0, method="bilinear")
        self.assertEqual(result, _node_value(1, 1))

    def test_interpolate_on_node_bicubic(self):
        result = interpolate_ntv2(self.grid, 1.0, -1.0, method="bicubic")
        self.assertEqual(result, _node_value(1, 1))

    def test_interpolate_outside_grid_returns_none_fields(self):
        result = interpolate_ntv2(self.grid, 50.0, -50.0, method="bicubic")
        self.assertEqual(result, (None, None, None, None))

    def test_interpolate_invalid_method_raises_value_error(self):
        with self.assertRaises(ValueError):
            interpolate_ntv2(self.grid, 1.0, -1.0, method="quadratic")


class TestInterpolateChoosesFinestOverlappingSubgrid(unittest.TestCase):
    """
    interpolate_ntv2 selects, among all subgrids whose extent contains the
    point, the one with the smallest lat_inc - it does not use the PARENT
    field at all. This builds a coarse 3x3 "parent" grid and a fine 4x4
    "child" grid with overlapping extents and distinct node values, and
    checks the fine grid's value wins.
    """

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmpdir.name, "synthetic_nested.gsb")

        coarse_nodes = [[(999.0, 999.0, 999.0, 999.0)] * 3 for _ in range(3)]
        coarse = Ntv2SubgridSpec(
            sub_name="COARSE",
            s_lat=0.0,
            n_lat=6 * 3600.0,
            e_long=0.0,
            w_long=6 * 3600.0,
            lat_inc=3 * 3600.0,
            long_inc=3 * 3600.0,
            nodes=coarse_nodes,
        )

        fine_nodes = [
            [(row * 10 + col + 0.5, -(row * 10 + col + 0.5), 1.0, 2.0) for col in range(4)]
            for row in range(4)
        ]
        fine = Ntv2SubgridSpec(
            sub_name="FINE",
            s_lat=0.0,
            n_lat=3 * 3600.0,
            e_long=0.0,
            w_long=3 * 3600.0,
            lat_inc=3600.0,
            long_inc=3600.0,
            nodes=fine_nodes,
            parent="COARSE",
        )

        write_ntv2_file(self.path, [coarse, fine])
        self.grid = read_ntv2_file(self.path)

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_both_subgrids_contain_the_point(self):
        self.assertEqual(len(self.grid.subgrids), 2)

    def test_finest_subgrid_value_is_used(self):
        result = interpolate_ntv2(self.grid, 1.0, -1.0, method="bicubic")
        self.assertEqual(result, (11.5, -11.5, 1.0, 2.0))


def _load_real_world_control_points():
    """
    Merges the committed, section-organised registry
    (ntv2_real_world_control_points.py) with an optional local overlay
    (ntv2_real_world_control_points_local.py, gitignored) for ad hoc grids
    that are only on your machine. Neither module is required to exist -
    a fresh checkout or CI simply gets zero control points, not an error.
    """
    points = []
    for module_name in (
        "geodepy.tests.ntv2_real_world_control_points",
        "geodepy.tests.ntv2_real_world_control_points_local",
    ):
        try:
            module = import_module(module_name)
        except ImportError:
            continue
        points.extend(getattr(module, "CONTROL_POINTS", []))
    return points


def _slug(text):
    return re.sub(r"\W+", "_", text).strip("_")


def _make_control_point_test(control_point):
    def test(self):
        resolved_path = resolve_grid_path(control_point.grid_path)
        if resolved_path is None:
            self.skipTest(
                f"set {GRIDS_DIR_ENV_VAR} to the folder containing "
                f"{control_point.grid_path} (see geodepy/tests/"
                f"ntv2_real_world_control_points.py) to run this check"
            )
        if not os.path.isfile(resolved_path):
            self.skipTest(f"grid file not found: {resolved_path}")

        grid = read_ntv2_file(resolved_path)
        result = interpolate_ntv2(
            grid, control_point.lat, control_point.lon, method=control_point.method
        )
        fields = dict(zip(("field_1", "field_2", "field_3", "field_4"), result))
        for field_name, expected_value in control_point.expected.items():
            with self.subTest(field=field_name):
                self.assertAlmostEqual(
                    fields[field_name], expected_value, places=control_point.places
                )

    test.__doc__ = f"{control_point.label} ({control_point.source})"
    return test


class TestNtv2RealWorldControlPoints(unittest.TestCase):
    """
    One generated test per entry in
    geodepy/tests/ntv2_real_world_control_points.py (committed, organised
    into sections per ICSM grid product - see
    https://github.com/icsm-au/transformation_grids) plus
    ntv2_real_world_control_points_local.py (optional, gitignored, for
    your own grids).

    These are opt-in by design: a bare filename entry only resolves to an
    actual file once you set GEODEPY_NTV2_GRIDS_DIR to wherever you've
    downloaded grids to, so a plain checkout/CI run sees every such entry
    skip rather than fail. A grid that's still missing even with the
    directory set is skipped individually too, so you only need to
    download the grids you actually want to check.
    """


_control_points = _load_real_world_control_points()
if _control_points:
    for _i, _cp in enumerate(_control_points):
        setattr(
            TestNtv2RealWorldControlPoints,
            f"test_{_i:02d}_{_slug(_cp.label)}",
            _make_control_point_test(_cp),
        )
else:

    def _test_no_control_points_registered(self):
        self.skipTest(
            "no entries in geodepy/tests/ntv2_real_world_control_points.py "
            "yet - see that file's docstring to add some"
        )

    TestNtv2RealWorldControlPoints.test_no_control_points_registered = (
        _test_no_control_points_registered
    )


if __name__ == "__main__":
    unittest.main()
