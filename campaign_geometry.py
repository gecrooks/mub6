"""Geometry-only tile-count and A100-hour estimates.

This module does not assert that a proposed half-width certifies. It turns
already measured rigorous widths into cover counts, keeping codimension and
anisotropy explicit so worst-axis cube pricing is easy to detect.
"""

import argparse
import json
import math


def axis_tile_count(length, half_width):
    if length < 0 or half_width <= 0:
        raise ValueError("length must be nonnegative and half-width positive")
    return max(1, math.ceil(length / (2.0 * half_width)))


def box_tile_count(lengths, half_widths, symmetry_factor=1.0):
    lengths = tuple(float(x) for x in lengths)
    half_widths = tuple(float(x) for x in half_widths)
    if len(lengths) != len(half_widths) or not lengths:
        raise ValueError("lengths and half_widths must have equal dimension")
    if symmetry_factor <= 0:
        raise ValueError("symmetry_factor must be positive")
    raw = math.prod(axis_tile_count(length, width)
                    for length, width in zip(lengths, half_widths))
    return math.ceil(raw / symmetry_factor)


def wall_band_comparison(tangent_lengths, band_width, normal_half_width,
                         tangent_half_width, seconds_per_tile=1.0,
                         symmetry_factor=1.0):
    """Compare anisotropic and worst-width cube covers of a wall band."""
    tangent_lengths = tuple(float(x) for x in tangent_lengths)
    lengths = tangent_lengths + (float(band_width),)
    anisotropic_h = (float(tangent_half_width),) * len(tangent_lengths) \
        + (float(normal_half_width),)
    isotropic_h = (float(normal_half_width),) * len(lengths)
    anisotropic = box_tile_count(lengths, anisotropic_h, symmetry_factor)
    isotropic = box_tile_count(lengths, isotropic_h, symmetry_factor)
    return {
        "lengths": list(lengths),
        "anisotropic_half_widths": list(anisotropic_h),
        "isotropic_half_widths": list(isotropic_h),
        "anisotropic_tiles": anisotropic,
        "isotropic_tiles": isotropic,
        "tile_reduction": isotropic / anisotropic,
        "anisotropic_a100_hours": anisotropic * seconds_per_tile / 3600.0,
        "isotropic_a100_hours": isotropic * seconds_per_tile / 3600.0,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tangent-lengths", default="1.5707963268,3.1415926536")
    parser.add_argument("--band-width", type=float, default=2.5e-5)
    parser.add_argument("--normal-h", type=float, default=1e-6)
    parser.add_argument("--tangent-h", type=float, default=2.5e-5)
    parser.add_argument("--seconds", type=float, default=6.0)
    parser.add_argument("--symmetry", type=float, default=1.0)
    args = parser.parse_args()
    lengths = tuple(float(x) for x in args.tangent_lengths.split(","))
    result = wall_band_comparison(
        lengths, args.band_width, args.normal_h, args.tangent_h,
        seconds_per_tile=args.seconds, symmetry_factor=args.symmetry,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
