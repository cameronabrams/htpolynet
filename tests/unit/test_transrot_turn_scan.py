"""

.. module:: test_transrot_turn_scan
   :synopsis: template placement turns the incoming piece to clear crowded sites

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import numpy as np
import pytest

from htpolynet.core.molecule import _min_distance, _rotate_about_axis


class TestRotateAboutAxis:
    def test_quarter_turn_about_z(self):
        got = _rotate_about_axis(np.array([[1.0, 0.0, 0.0]]), np.zeros(3), np.array([0.0, 0.0, 1.0]), np.pi / 2)
        assert np.allclose(got, [[0.0, 1.0, 0.0]])

    def test_points_on_the_axis_stay_put(self):
        pts = np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 5.0]])
        got = _rotate_about_axis(pts, np.array([1.0, 2.0, 0.0]), np.array([0.0, 0.0, 2.0]), 1.234)
        assert np.allclose(got, pts)

    def test_axis_through_an_offset_point(self):
        got = _rotate_about_axis(np.array([[2.0, 0.0, 0.0]]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0]), np.pi)
        assert np.allclose(got, [[0.0, 0.0, 0.0]])

    def test_distances_are_preserved(self):
        rng = np.random.default_rng(1)
        pts = rng.normal(size=(6, 3))
        got = _rotate_about_axis(pts, rng.normal(size=3), rng.normal(size=3), 0.7)
        d0 = np.linalg.norm(pts[:, None] - pts[None], axis=2)
        d1 = np.linalg.norm(got[:, None] - got[None], axis=2)
        assert np.allclose(d0, d1)

    def test_zero_turn_is_a_copy(self):
        pts = np.ones((2, 3))
        got = _rotate_about_axis(pts, np.zeros(3), np.array([1.0, 0.0, 0.0]), 0.0)
        got[0, 0] = 9.0
        assert pts[0, 0] == 1.0


class TestMinDistance:
    def test_closest_pair(self):
        a = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
        b = np.array([[0.0, 3.0, 0.0], [5.0, 0.5, 0.0]])
        assert _min_distance(a, b) == pytest.approx(0.5)

    def test_empty_is_infinite(self):
        assert _min_distance(np.empty((0, 3)), np.zeros((1, 3))) == np.inf


class TestTurnScanClearsAClash:
    def test_the_best_turn_moves_a_substituent_off_an_occupied_site(self):
        # one arm of a piece bonded along the z axis; a fixed atom sits where the arm points at 0 deg
        piece = np.array([[0.15, 0.0, 0.35]])
        fixed = np.array([[0.15, 0.0, 0.36]])
        scores = []
        for theta in np.deg2rad(np.arange(0, 360, 10)):
            scores.append((_min_distance(_rotate_about_axis(piece, np.zeros(3), np.array([0, 0, 1.0]), theta), fixed),
                           theta))
        best = max(scores)
        assert scores[0][0] < 0.02
        assert best[0] > 0.25
        assert np.isclose(best[1], np.pi)
