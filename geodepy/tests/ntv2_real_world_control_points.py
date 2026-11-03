"""
Committed registry of independently verified NTv2 control points,
consumed by TestNtv2RealWorldControlPoints in test_ntv2reader.py.

Organised into one section per grid published in ICSM's
transformation_grids repository (https://github.com/icsm-au/transformation_grids),
which the README there describes as "the most commonly used transformation
grids in Australia":

    - AGD66 -> GDA94            A66_National_13_09_01.gsb
    - AGD84 -> GDA94            National_84_02_07_01.gsb
    - GDA94 -> GDA2020 conformal             GDA94_GDA2020_conformal.gsb
    - GDA94 -> GDA2020 conformal+distortion  GDA94_GDA2020_conformal_and_distortion.gsb

This file only contains control point *values* (lat/lon and their
independently-sourced expected shifts), not the grids themselves (tens to
hundreds of MB, licensed separately under CC BY 4.0 by ICSM) - it's safe
to commit. Each Ntv2ControlPoint.grid_path here is a bare filename, which
geodepy.tests.ntv2_builder.resolve_grid_path() looks up under the
GEODEPY_NTV2_GRIDS_DIR environment variable at test time. That makes
these tests opt-in: with the variable unset, every entry here simply
skips (safe for a default checkout or CI); set it once you've downloaded
the grids you want to check against, e.g. (PowerShell):

    $env:GEODEPY_NTV2_GRIDS_DIR = "D:\\grids\\icsm_transformation_grids"
    python -m unittest geodepy.tests.test_ntv2reader.TestNtv2RealWorldControlPoints -v

A grid you haven't downloaded is skipped individually - you don't need
all four to run some of these.

For your own, non-ICSM grids (e.g. an AUSGeoid release) kept on an
absolute, machine-specific path, add them to the optional, gitignored
geodepy/tests/ntv2_real_world_control_points_local.py instead, rather
than to this file.

IMPORTANT: every `expected` value below must come from a source
independent of geodepy.ntv2reader itself (a technical manual, an official
online calculator, a published test point) - never from running
interpolate_ntv2() and copying its output, which can't catch a regression
in the function under test. Record that source in `source`.
"""

from geodepy.tests.ntv2_builder import Ntv2ControlPoint

CONTROL_POINTS = []

# ---------------------------------------------------------------------------
# AGD66 -> GDA94  |  A66_National_13_09_01.gsb  |  EPSG 1803
# AGD66 (EPSG 4202) <-> GDA94 (EPSG 4283), full national coverage.
# https://github.com/icsm-au/transformation_grids
# ---------------------------------------------------------------------------
CONTROL_POINTS += [
    # Ntv2ControlPoint(
    #     grid_path="A66_National_13_09_01.gsb",
    #     lat=-35.3081,
    #     lon=149.1244,
    #     expected={"field_1": 0.0, "field_2": 0.0},  # TODO: source an AGD66<->GDA94 test point
    #     source="TODO: e.g. a published AGD66/GDA94 control point pair",
    #     label="A66_National_example",
    # ),
]

# ---------------------------------------------------------------------------
# AGD84 -> GDA94  |  National_84_02_07_01.gsb  |  EPSG 1804
# AGD84 (EPSG 4203) <-> GDA94 (EPSG 4283), WA/SA/QLD coverage.
# https://github.com/icsm-au/transformation_grids
# ---------------------------------------------------------------------------
CONTROL_POINTS += [
    # Ntv2ControlPoint(
    #     grid_path="National_84_02_07_01.gsb",
    #     lat=-31.9523,
    #     lon=115.8613,
    #     expected={"field_1": 0.0, "field_2": 0.0},  # TODO: source an AGD84<->GDA94 test point
    #     source="TODO: e.g. a published AGD84/GDA94 control point pair",
    #     label="National_84_example",
    # ),
]

# ---------------------------------------------------------------------------
# GDA94 -> GDA2020 (conformal)  |  GDA94_GDA2020_conformal.gsb  |  EPSG 8446
# GDA94 (EPSG 4283) <-> GDA2020 (EPSG 7844), full national coverage.
# A pure 7-parameter conformal (Helmert) transformation, gridded - mainly
# plate-tectonic motion, ICSM's README quotes ~1.8 m NNE.
# https://github.com/icsm-au/transformation_grids
# ---------------------------------------------------------------------------
CONTROL_POINTS += [
    Ntv2ControlPoint(
        grid_path="GDA94_GDA2020_conformal.gsb",
        lat=-23.670123894062947,
        lon=133.8855132902714,
        # ALIC (Alice Springs). Derived from the GDA94/GDA2020 Cartesian
        # coordinates already cited in geodepy/tests/test_transform.py
        # (GDA2020 Technical Manual v1.2 Sect 3.1.1), converted to
        # geographic via geodepy.convert.xyz2llh - not via interpolate_ntv2.
        # shift_lat = (lat_gda2020 - lat_gda94) * 3600
        # shift_lon = (lon_gda94 - lon_gda2020) * 3600
        # Resulting ~1.75 m NNE shift matches the ICSM README's ~1.8 m NNE
        # figure for this grid, as an extra sanity check on the derivation.
        expected={"field_1": 0.04952, "field_2": -0.029947},
        source=(
            "GDA2020 Technical Manual v1.2 Sect 3.1.1 ALIC GDA94/GDA2020 "
            "Cartesian coordinates (as used in test_transform.py), "
            "converted to geographic via geodepy.convert.xyz2llh"
        ),
        label="GDA94_GDA2020_conformal_ALIC",
        places=3,  # ~1 mm; pure-conformal grid, minimal discretisation error expected
    ),
]

# ---------------------------------------------------------------------------
# GDA94 -> GDA2020 (conformal + distortion)
# GDA94_GDA2020_conformal_and_distortion.gsb  |  EPSG 8447
# GDA94 (EPSG 4283) <-> GDA2020 (EPSG 7844), full national coverage.
# Conformal transformation plus a regional distortion component
# (reference frame realisation + ground movement, up to ~0.5 m) - values
# here will generally NOT match the pure-conformal grid above.
# https://github.com/icsm-au/transformation_grids
# ---------------------------------------------------------------------------
CONTROL_POINTS += [
    # Ntv2ControlPoint(
    #     grid_path="GDA94_GDA2020_conformal_and_distortion.gsb",
    #     lat=-23.670123894062947,
    #     lon=133.8855132902714,
    #     expected={"field_1": 0.0, "field_2": 0.0},  # TODO: source a distortion-grid test point
    #     source="TODO: e.g. GDA2020 Technical Manual Table 3.4 or GA's online transformation tool",
    #     label="GDA94_GDA2020_conformal_and_distortion_ALIC",
    # ),
]
