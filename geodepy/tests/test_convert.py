import os
import unittest
import datetime
import numpy as np
from geodepy.fileio import read_dnacoord
from geodepy.constants import grs80, ans, isg
from geodepy.convert import (
    hp2dec,
    DMSAngle,
    DDMAngle,
    yyyydoy_to_date,
    date_to_yyyydoy,
    grid2geo,
    hp2dec_v,
    geo2grid,
    llh2xyz,
    xyz2llh,
    xyz2grid,
    grid2xyz,
)


class TestConvert(unittest.TestCase):
    def test_date_to_yyyydoy(self):
        self.assertEqual(date_to_yyyydoy(datetime.date(2020, 1, 4)), "2020.004")
        self.assertEqual(date_to_yyyydoy(datetime.date(2020, 10, 12)), "2020.286")
        self.assertEqual(date_to_yyyydoy(datetime.date(1998, 4, 7)), "1998.097")
        self.assertEqual(date_to_yyyydoy(datetime.date(2000, 11, 22)), "2000.327")
        self.assertEqual(date_to_yyyydoy(datetime.date(2008, 2, 29)), "2008.060")
        with self.assertRaises(AttributeError):
            date_to_yyyydoy("a")
        with self.assertRaises(AttributeError):
            date_to_yyyydoy("2020123")

    def test_yyyydoy_to_date(self):
        self.assertEqual(yyyydoy_to_date("2020.004"), datetime.date(2020, 1, 4))
        self.assertEqual(yyyydoy_to_date("2020.286"), datetime.date(2020, 10, 12))
        self.assertEqual(yyyydoy_to_date("1998.097"), datetime.date(1998, 4, 7))
        self.assertEqual(yyyydoy_to_date("2000.327"), datetime.date(2000, 11, 22))
        self.assertEqual(yyyydoy_to_date("2008.060"), datetime.date(2008, 2, 29))
        self.assertEqual(yyyydoy_to_date("2020004"), datetime.date(2020, 1, 4))
        self.assertEqual(yyyydoy_to_date("2020286"), datetime.date(2020, 10, 12))
        self.assertEqual(yyyydoy_to_date("1998097"), datetime.date(1998, 4, 7))
        self.assertEqual(yyyydoy_to_date("2000327"), datetime.date(2000, 11, 22))
        self.assertEqual(yyyydoy_to_date("2008060"), datetime.date(2008, 2, 29))
        with self.assertRaises(ValueError):
            yyyydoy_to_date("a")
        with self.assertRaises(ValueError):
            yyyydoy_to_date("20201234")
        with self.assertRaises(ValueError):
            yyyydoy_to_date("2020.1234")
        with self.assertRaises(ValueError):
            yyyydoy_to_date("202012")
        with self.assertRaises(ValueError):
            yyyydoy_to_date("2020.12")

    def test_geo_grid_transform_interoperability(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))

        test_geo_coords = np.genfromtxt(
            os.path.join(abs_path, "resources/Test_Conversion_Geo.csv"),
            delimiter=",",
            dtype="S4,f8,f8",
            names=["site", "lat", "lon"],
        )

        test_grid_coords = np.genfromtxt(
            os.path.join(abs_path, "resources/Test_Conversion_Grid.csv"),
            delimiter=",",
            dtype="S4,i4,f8,f8",
            names=["site", "zone", "east", "north"],
        )

        geoed_grid = np.array(
            [grid2geo(*x) for x in test_grid_coords[["zone", "east", "north"]]]
        )
        np.testing.assert_almost_equal(
            geoed_grid[:, :2],
            hp2dec_v(np.array(test_geo_coords[["lat", "lon"]].tolist())),
            decimal=8,
        )

        gridded_geo = np.stack(
            [
                geo2grid(*x)
                for x in hp2dec_v(np.array(test_geo_coords[["lat", "lon"]].tolist()))
            ]
        )
        np.testing.assert_almost_equal(
            gridded_geo[:, 2:4].astype(float),
            np.array(test_grid_coords[["east", "north"]].tolist()),
            decimal=3,
        )

    def test_geo_grid_transform_interoperability_isg(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))
        test_geo_coords = np.genfromtxt(
            os.path.join(abs_path, "resources/Test_Conversion_ISG_Geo.csv"),
            delimiter=",",
            dtype="S4,f8,f8",
            names=["site", "lat", "lon"],
        )

        test_grid_coords = np.genfromtxt(
            os.path.join(abs_path, "resources/Test_Conversion_ISG_Grid.csv"),
            delimiter=",",
            dtype="S4,i4,f8,f8",
            names=["site", "zone", "east", "north"],
        )

        geoed_grid = np.array(
            [
                grid2geo(*x, ellipsoid=ans, prj=isg)
                for x in test_grid_coords[["zone", "east", "north"]]
            ]
        )
        np.testing.assert_almost_equal(
            geoed_grid[:, :2],
            np.array(test_geo_coords[["lat", "lon"]].tolist()),
            decimal=8,
        )

        gridded_geo = np.stack(
            [
                geo2grid(*x, ellipsoid=ans, prj=isg)
                for x in np.array(test_geo_coords[["lat", "lon"]].tolist())
            ]
        )
        np.testing.assert_almost_equal(
            gridded_geo[:, 2:4].astype(float),
            np.array(test_grid_coords[["east", "north"]].tolist()),
            decimal=3,
        )

    def test_llh2xyz(self):

        # Test of single point
        x, y, z = llh2xyz(hp2dec(-37.482667598), hp2dec(144.581644114), 39.6514)
        self.assertAlmostEqual(x, -4131654.2815, 3)
        self.assertAlmostEqual(y, 2896107.9738, 3)
        self.assertAlmostEqual(z, -3888601.3067, 3)

        # Test DMSAngle input
        x, y, z = llh2xyz(
            DMSAngle(-37, 48, 26.67598), DMSAngle(144, 58, 16.44114), 39.6514
        )
        self.assertAlmostEqual(x, -4131654.2815, 3)
        self.assertAlmostEqual(y, 2896107.9738, 3)
        self.assertAlmostEqual(z, -3888601.3067, 3)

        # Test DDMAngle input
        x, y, z = llh2xyz(DDMAngle(-37, 48.4445996), DDMAngle(144, 58.274019), 39.6514)
        self.assertAlmostEqual(x, -4131654.2815, 3)
        self.assertAlmostEqual(y, 2896107.9738, 3)
        self.assertAlmostEqual(z, -3888601.3067, 3)

        # Tests comparing results from DynAdjust Adjusted Coordinate File
        # including 109 Points Across Australia
        abs_path = os.path.abspath(os.path.dirname(__file__))

        testdata = read_dnacoord(
            os.path.join(abs_path, "resources/natadjust_rvs_example.dat")
        )
        for coord in testdata:
            coord.converthptodd()
            xcomp, ycomp, zcomp = llh2xyz(coord.lat, coord.long, coord.ell_ht)
            assert abs(coord.x - xcomp) < 2e-4
            assert abs(coord.y - ycomp) < 2e-4
            assert abs(coord.z - zcomp) < 2e-4

        # Test point on equator uses the supplied ellipsoid
        x, y, z = llh2xyz(0, 0, 0, ellipsoid=ans)
        self.assertAlmostEqual(x, ans.semimaj, 6)
        self.assertAlmostEqual(y, 0, 6)
        self.assertAlmostEqual(z, 0, 6)

    def test_xyz2llh(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))

        testdata = read_dnacoord(
            os.path.join(abs_path, "resources/natadjust_rvs_example.dat")
        )
        for coord in testdata:
            coord.converthptodd()
            latcomp, longcomp, ell_htcomp = xyz2llh(coord.x, coord.y, coord.z)
            assert abs(latcomp - coord.lat) < 2e-9
            assert abs(longcomp - coord.long) < 2e-9
            assert abs(ell_htcomp - coord.ell_ht) < 1e-4

    def test_geo2grid(self):
        # Single Point Test
        hem, zone, east, north, psf, grid_conv = geo2grid(
            hp2dec(-37.482667598), hp2dec(144.581644114)
        )
        self.assertEqual(hem, "South")
        self.assertEqual(zone, 55)
        self.assertAlmostEqual(east, 321405.5592, 3)
        self.assertAlmostEqual(north, 5813614.1613, 3)
        self.assertAlmostEqual(psf, 0.99999287, 8)
        self.assertAlmostEqual(grid_conv, -1.2439811331, 9)

        # Test DMSAngle Input
        (hem, zone, east, north, psf, grid_conv) = geo2grid(
            DMSAngle(-37, 48, 26.67598), DMSAngle(144, 58, 16.44114)
        )
        self.assertEqual(hem, "South")
        self.assertEqual(zone, 55)
        self.assertAlmostEqual(east, 321405.5592, 3)
        self.assertAlmostEqual(north, 5813614.1613, 3)
        self.assertAlmostEqual(psf, 0.99999287, 8)
        self.assertAlmostEqual(grid_conv, -1.2439811331, 9)

        # Test DDMAngle Input
        (hem, zone, east, north, psf, grid_conv) = geo2grid(
            DDMAngle(-37, 48.4445997), DDMAngle(144, 58.274019)
        )
        self.assertEqual(hem, "South")
        self.assertEqual(zone, 55)
        self.assertAlmostEqual(east, 321405.5592, 3)
        self.assertAlmostEqual(north, 5813614.1613, 3)
        self.assertAlmostEqual(psf, 0.99999287, 8)
        self.assertAlmostEqual(grid_conv, -1.2439811331, 9)

        abs_path = os.path.abspath(os.path.dirname(__file__))

        # Test various coordinates in Australia
        testdata = read_dnacoord(
            os.path.join(abs_path, "resources/natadjust_rvs_example.dat")
        )
        for coord in testdata:
            coord.converthptodd()
            hem, zonecomp, eastcomp, northcomp, psf, grid_conv = geo2grid(
                coord.lat, coord.long
            )
            self.assertEqual(zonecomp, coord.zone)
            self.assertLess(abs(eastcomp - coord.easting), 4e-4)
            self.assertLess((northcomp - coord.northing), 4e-4)

        # Test North and South Hemisphere Output
        north_ex = (DMSAngle(34, 57, 00.79653).dec(), DMSAngle(117, 48, 36.68783).dec())
        south_ex = (
            DMSAngle(-34, 57, 00.79653).dec(),
            DMSAngle(117, 48, 36.68783).dec(),
        )
        north_grid = geo2grid(north_ex[0], north_ex[1])
        south_grid = geo2grid(south_ex[0], south_ex[1])
        self.assertEqual(north_grid[0], "North")
        self.assertEqual(south_grid[0], "South")
        self.assertEqual(north_grid[1], south_grid[1])  # Zone
        self.assertEqual(north_grid[2], south_grid[2])  # Easting
        self.assertEqual(north_grid[3], 10000000 - south_grid[3])  # Northing
        self.assertEqual(north_grid[4], south_grid[4])  # PSF
        self.assertEqual(north_grid[5], -south_grid[5])  # Grid Convergence

        # Test Input Validation
        with self.assertRaises(ValueError):
            geo2grid(0, 45, -1)
        with self.assertRaises(ValueError):
            geo2grid(0, 45, 61)
        with self.assertRaises(ValueError):
            geo2grid(-81, 45, 0)
        with self.assertRaises(ValueError):
            geo2grid(85, 45, 0)
        with self.assertRaises(ValueError):
            geo2grid(0, -181, 0)
        with self.assertRaises(ValueError):
            geo2grid(0, 181, 0)
        # test ValueError raised if not a valid ISG zone
        with self.assertRaises(ValueError):
            grid2geo(zone=530, east=300000, north=1500000, ellipsoid=ans, prj=isg)
        # test UserWarning raised if prj/ellipsoid combination isn't recommended
        with self.assertWarns(UserWarning):
            grid2geo(zone=551, east=300000, north=1500000, ellipsoid=grs80, prj=isg)

    def test_grid2geo(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))

        testdata = read_dnacoord(
            os.path.join(abs_path, "resources/natadjust_rvs_example.dat")
        )
        for coord in testdata:
            coord.converthptodd()
            latcomp, longcomp, psf, grid_conv = grid2geo(
                coord.zone, coord.easting, coord.northing
            )
            self.assertLess(abs(latcomp - coord.lat), 5e-9)
            self.assertLess(abs(longcomp - coord.long), 5e-9)

        # Test North and South Hemisphere Output
        north_ex = (50, 573976.8747, 3867822.4539, "North")
        south_ex = (50, 573976.8747, 6132177.5461, "South")
        north_geo = grid2geo(*north_ex)
        south_geo = grid2geo(*south_ex)
        self.assertEqual(north_geo[0], -south_geo[0])
        self.assertEqual(north_geo[1], south_geo[1])
        self.assertEqual(north_geo[2], south_geo[2])
        self.assertEqual(north_geo[3], -south_geo[3])

        # Test Input Validation
        with self.assertRaises(ValueError):
            grid2geo(-1, 0, 500000)
        with self.assertRaises(ValueError):
            grid2geo(61, 0, 500000)
        with self.assertRaises(ValueError):
            grid2geo(0, -2830001, 500000)
        with self.assertRaises(ValueError):
            grid2geo(0, 3830001, 500000)
        with self.assertRaises(ValueError):
            grid2geo(0, 0, -1)
        with self.assertRaises(ValueError):
            grid2geo(0, 0, 10000001)
        with self.assertRaises(ValueError):
            grid2geo(0, 0, 500000, "fail")
        # test ValueError raised if not a valid ISG zone
        with self.assertRaises(ValueError):
            grid2geo(zone=530, east=300000, north=1500000, ellipsoid=ans, prj=isg)
        # test UserWarning raised if prj/ellipsoid combination isn't recommended
        with self.assertWarns(UserWarning):
            grid2geo(zone=551, east=300000, north=1500000, ellipsoid=grs80, prj=isg)

    def test_xyz2grid(self):
        # Single Point Test
        hem, zone, east, north, ellht, psf, grid_conv = xyz2grid(
            -4131654.2815, 2896107.9738, -3888601.3067
        )
        self.assertEqual(hem, "South")
        self.assertEqual(zone, 55)
        self.assertAlmostEqual(east, 321405.5592, 3)
        self.assertAlmostEqual(north, 5813614.1613, 3)
        self.assertAlmostEqual(ellht, 39.6514, 3)
        self.assertAlmostEqual(psf, 0.99999287, 8)
        self.assertAlmostEqual(grid_conv, -1.2439811331, 9)

        # Tests comparing results from DynAdjust Adjusted Coordinate File
        # including 109 Points Across Australia
        abs_path = os.path.abspath(os.path.dirname(__file__))
        testdata = read_dnacoord(
            os.path.join(abs_path, "resources/natadjust_rvs_example.dat")
        )
        for coord in testdata:
            coord.converthptodd()
            hem, zonecomp, eastcomp, northcomp, ellhtcomp, psf, grid_conv = xyz2grid(
                coord.x, coord.y, coord.z
            )
            self.assertEqual(zonecomp, coord.zone)
            self.assertLess(abs(eastcomp - coord.easting), 4e-4)
            self.assertLess(abs(northcomp - coord.northing), 4e-4)
            self.assertLess(abs(ellhtcomp - coord.ell_ht), 2e-4)

        # Test North Hemisphere Output
        x, y, z = llh2xyz(DMSAngle(34, 57, 00.79653), DMSAngle(117, 48, 36.68783))
        north_grid = xyz2grid(x, y, z)
        self.assertEqual(north_grid[0], "North")
        self.assertEqual(north_grid[1], 50)
        self.assertAlmostEqual(north_grid[2], 573976.8747, 3)
        self.assertAlmostEqual(north_grid[3], 3867822.4539, 3)

        # Test ISG projection with ANS ellipsoid
        x, y, z = llh2xyz(-33.8, 151.2, 50, ellipsoid=ans)
        isg_grid = xyz2grid(x, y, z, ellipsoid=ans, prj=isg)
        geo_grid = geo2grid(-33.8, 151.2, ellipsoid=ans, prj=isg)
        self.assertEqual(isg_grid[1], geo_grid[1])  # Zone
        self.assertAlmostEqual(isg_grid[2], geo_grid[2], 3)  # Easting
        self.assertAlmostEqual(isg_grid[3], geo_grid[3], 3)  # Northing
        self.assertAlmostEqual(isg_grid[4], 50, 3)  # Ellipsoid Height

    def test_grid2xyz(self):
        # Single Point Test
        x, y, z = grid2xyz(55, 321405.5592, 5813614.1613, 39.6514)
        self.assertAlmostEqual(x, -4131654.2815, 3)
        self.assertAlmostEqual(y, 2896107.9738, 3)
        self.assertAlmostEqual(z, -3888601.3067, 3)

        # Tests comparing results from DynAdjust Adjusted Coordinate File
        # including 109 Points Across Australia
        abs_path = os.path.abspath(os.path.dirname(__file__))
        testdata = read_dnacoord(
            os.path.join(abs_path, "resources/natadjust_rvs_example.dat")
        )
        for coord in testdata:
            coord.converthptodd()
            xcomp, ycomp, zcomp = grid2xyz(
                coord.zone, coord.easting, coord.northing, coord.ell_ht
            )
            self.assertLess(abs(xcomp - coord.x), 3e-4)
            self.assertLess(abs(ycomp - coord.y), 3e-4)
            self.assertLess(abs(zcomp - coord.z), 3e-4)

        # Test North Hemisphere Input
        x, y, z = grid2xyz(50, 573976.8747, 3867822.4539, hemisphere="North")
        xref, yref, zref = llh2xyz(
            DMSAngle(34, 57, 00.79653), DMSAngle(117, 48, 36.68783)
        )
        self.assertAlmostEqual(x, xref, 3)
        self.assertAlmostEqual(y, yref, 3)
        self.assertAlmostEqual(z, zref, 3)

        # Test round trip through xyz2grid
        xyz = (-4131654.2815, 2896107.9738, -3888601.3067)
        hem, zone, east, north, ellht, psf, grid_conv = xyz2grid(*xyz)
        np.testing.assert_allclose(
            grid2xyz(zone, east, north, ellht, hemisphere=hem), xyz, atol=1e-3
        )

        # Test round trip with ISG projection and ANS ellipsoid
        xyz = llh2xyz(-33.8, 151.2, 50, ellipsoid=ans)
        hem, zone, east, north, ellht, psf, grid_conv = xyz2grid(
            *xyz, ellipsoid=ans, prj=isg
        )
        np.testing.assert_allclose(
            grid2xyz(
                zone, east, north, ellht, hemisphere=hem, ellipsoid=ans, prj=isg
            ),
            xyz,
            atol=1e-3,
        )

        # Test hemisphere, ellipsoid and prj are keyword-only
        with self.assertRaises(TypeError):
            grid2xyz(50, 573976.8747, 3867822.4539, 0, "North")

        # Test Input Validation passes through from grid2geo
        with self.assertRaises(ValueError):
            grid2xyz(61, 500000, 5000000)
        with self.assertRaises(ValueError):
            grid2xyz(55, 500000, 5000000, hemisphere="fail")
        with self.assertRaises(ValueError):
            grid2xyz(530, 300000, 1500000, ellipsoid=ans, prj=isg)


if __name__ == "__main__":
    unittest.main()
