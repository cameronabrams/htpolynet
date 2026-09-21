"""Geometric helpers for putting an atom or a group somewhere it fits.

Used when a bond has to be formed between atoms that no existing geometry places
next to each other: postcure repair relocating a -C#N cap, and the template build
positioning one reactant against another for an addition, where there are no
sacrificial hydrogens to align on.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging

import numpy as np

logger = logging.getLogger(__name__)


def sphere_directions(n=48):
    """``n`` roughly uniform unit vectors on the sphere, deterministically.

    A Fibonacci spiral: no randomness, so a rebuild of the same system makes the
    same choices.

    Args:
        n (int): how many directions to generate

    Returns:
        numpy.ndarray: (n, 3) array of unit vectors
    """
    k = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * k / n)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * k
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


def rotate_about_axis(xyz, point, axis, theta):
    """Rotates points by theta about the line through point along axis (Rodrigues).

    Args:
        xyz (numpy.ndarray): (n, 3) positions
        point (numpy.ndarray): a point on the axis
        axis (numpy.ndarray): direction of the axis, normalized here
        theta (float): angle in radians

    Returns:
        numpy.ndarray: rotated (n, 3) positions
    """
    if theta == 0.0:
        return xyz.copy()
    k = np.asarray(axis, dtype=float)
    k = k / np.linalg.norm(k)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    R = np.eye(3) + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)
    return (np.asarray(xyz, dtype=float) - point) @ R.T + point


def min_distance(a, b):
    """Smallest distance between any point of a and any point of b; inf if either is empty."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return np.inf
    diff = a[:, np.newaxis, :] - b[np.newaxis, :, :]
    return float(np.sqrt((diff * diff).sum(axis=2)).min())


def place_group(group, attach, anchor, fixed, bond_length, n_directions=48, n_turns=24):
    """Positions a rigid group so one of its atoms sits a bond length from an anchor atom.

    Two searches, both coarse and deterministic: where to put the attaching atom --
    a direction on the sphere around the anchor -- and how to turn the group about
    that direction once placed.  Each candidate is scored by its closest approach to
    the atoms that stay put, and the best is returned.

    Args:
        group (numpy.ndarray): (n, 3) positions of the group to move
        attach (int): index within group of the atom that will carry the new bond
        anchor (array-like): position of the atom the new bond reaches from
        fixed (numpy.ndarray): (m, 3) positions that must be cleared
        bond_length (float): distance to leave between the attaching atom and the anchor
        n_directions (int): sphere directions to try
        n_turns (int): turns about the chosen direction to try

    Returns:
        tuple: (positions, clearance) -- the placed group, and its closest approach to
            ``fixed`` excluding the attaching atom itself
    """
    group = np.asarray(group, dtype=float)
    anchor = np.asarray(anchor, dtype=float)
    others = [i for i in range(len(group)) if i != attach]
    best = None
    for u in sphere_directions(n_directions):
        target = anchor + bond_length * u
        shifted = group + (target - group[attach])
        for theta in np.linspace(0.0, 2.0 * np.pi, n_turns, endpoint=False):
            trial = rotate_about_axis(shifted, target, u, theta)
            clearance = min_distance(trial[others], fixed) if others else np.inf
            if best is None or clearance > best[1]:
                best = (trial, clearance)
    return best
