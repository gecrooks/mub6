"""Canonical campaign manifest suitable for an independent rechecker."""

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

from certificate_result import CertificateGrade
from ledger_bits import bits_float, float_bits


SCHEMA = "mub6-campaign-manifest-v2-cells"


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class LedgerShard:
    name: str
    sha256: str

    def __post_init__(self):
        if not self.name or Path(self.name).is_absolute() or ".." in Path(self.name).parts:
            raise ValueError("ledger shard name must be a safe relative path")
        if len(self.sha256) != 64 or any(c not in "0123456789abcdef"
                                        for c in self.sha256):
            raise ValueError("ledger shard digest must be lowercase SHA-256")

    def as_dict(self):
        return {"name": self.name, "sha256": self.sha256}


@dataclass(frozen=True)
class TransverseCell:
    line: tuple[float, float]
    phi_bounds: tuple[float, float]
    lambda_bounds: tuple[float, float]

    def __post_init__(self):
        if len(self.line) != 2 or len(self.phi_bounds) != 2 \
                or len(self.lambda_bounds) != 2:
            raise ValueError("transverse cell fields must be pairs")
        if self.phi_bounds[0] > self.phi_bounds[1] \
                or self.lambda_bounds[0] > self.lambda_bounds[1]:
            raise ValueError("transverse cell bounds are reversed")
        if not (self.phi_bounds[0] <= self.line[0] <= self.phi_bounds[1]) \
                or not (self.lambda_bounds[0] <= self.line[1]
                        <= self.lambda_bounds[1]):
            raise ValueError("transverse line lies outside its cell")

    def as_dict(self):
        return {
            "line_bits": [float_bits(value) for value in self.line],
            "phi_bounds_bits": [float_bits(value)
                                for value in self.phi_bounds],
            "lambda_bounds_bits": [float_bits(value)
                                   for value in self.lambda_bounds],
        }

    @classmethod
    def from_dict(cls, value):
        return cls(tuple(bits_float(x) for x in value["line_bits"]),
                   tuple(bits_float(x) for x in value["phi_bounds_bits"]),
                   tuple(bits_float(x)
                         for x in value["lambda_bounds_bits"]))


@dataclass(frozen=True)
class CampaignManifest:
    domain: tuple[tuple[float, float], ...]
    symmetry: str
    symmetry_factor: int
    required_grade: CertificateGrade
    ledger_shards: tuple[LedgerShard, ...]
    transverse_cells: tuple[TransverseCell, ...] = ()

    def __post_init__(self):
        if len(self.domain) != 3 or any(len(axis) != 2 for axis in self.domain):
            raise ValueError("campaign domain must have three intervals")
        if any(lo > hi for lo, hi in self.domain):
            raise ValueError("campaign domain interval is reversed")
        if not self.symmetry.strip() or self.symmetry_factor < 1:
            raise ValueError("campaign symmetry metadata is invalid")
        object.__setattr__(self, "required_grade",
                           CertificateGrade(self.required_grade))
        object.__setattr__(self, "ledger_shards", tuple(self.ledger_shards))
        names = [shard.name for shard in self.ledger_shards]
        if len(names) != len(set(names)):
            raise ValueError("campaign ledger shard names must be unique")
        cells = tuple(self.transverse_cells)
        object.__setattr__(self, "transverse_cells", cells)
        lines = [cell.line for cell in cells]
        if len(lines) != len(set(lines)):
            raise ValueError("transverse cell lines must be unique")
        if cells:
            phi_parts = sorted(set(cell.phi_bounds for cell in cells))
            lam_parts = sorted(set(cell.lambda_bounds for cell in cells))
            self._check_axis_partition(phi_parts, self.domain[1], "phi")
            self._check_axis_partition(lam_parts, self.domain[2], "lambda")
            actual = {(cell.phi_bounds, cell.lambda_bounds) for cell in cells}
            expected = {(phi, lam) for phi in phi_parts for lam in lam_parts}
            if actual != expected or len(cells) != len(expected):
                raise ValueError("transverse cells must form a Cartesian grid")

    @staticmethod
    def _check_axis_partition(parts, domain, name):
        if not parts or parts[0][0] != domain[0] \
                or parts[-1][1] != domain[1]:
            raise ValueError(f"{name} cells do not reach domain boundary")
        if any(left[1] != right[0] for left, right in zip(parts, parts[1:])):
            raise ValueError(f"{name} cells have a gap or overlap")

    def _unsigned_dict(self):
        return {
            "schema": SCHEMA,
            "domain_bits": [[float_bits(lo), float_bits(hi)]
                            for lo, hi in self.domain],
            "symmetry": self.symmetry,
            "symmetry_factor": self.symmetry_factor,
            "required_grade": self.required_grade.name,
            "ledger_shards": [shard.as_dict() for shard in self.ledger_shards],
            "transverse_cells": [cell.as_dict()
                                 for cell in self.transverse_cells],
        }

    @property
    def manifest_id(self):
        encoded = json.dumps(self._unsigned_dict(), sort_keys=True,
                             separators=(",", ":"), allow_nan=False)
        return hashlib.sha256(encoded.encode("utf-8")).hexdigest()

    def as_dict(self):
        return {**self._unsigned_dict(), "manifest_id": self.manifest_id}

    @classmethod
    def from_dict(cls, value):
        if value.get("schema") != SCHEMA:
            raise ValueError("unknown campaign manifest schema")
        manifest = cls(
            domain=tuple(tuple(bits_float(x) for x in axis)
                         for axis in value["domain_bits"]),
            symmetry=value["symmetry"],
            symmetry_factor=int(value["symmetry_factor"]),
            required_grade=CertificateGrade[value["required_grade"]],
            ledger_shards=tuple(LedgerShard(item["name"], item["sha256"])
                                for item in value["ledger_shards"]),
            transverse_cells=tuple(TransverseCell.from_dict(item)
                                   for item in value.get("transverse_cells", ())),
        )
        if value.get("manifest_id") != manifest.manifest_id:
            raise ValueError("campaign manifest digest mismatch")
        return manifest

    def verify_shards(self, base_directory):
        base = Path(base_directory)
        failures = []
        for shard in self.ledger_shards:
            path = base / shard.name
            try:
                actual = file_sha256(path)
            except OSError as error:
                failures.append(f"{shard.name}: {error}")
                continue
            if actual != shard.sha256:
                failures.append(f"{shard.name}: digest mismatch")
        return tuple(failures)


def build_manifest(domain, symmetry, symmetry_factor, required_grade,
                   shard_paths, *, base_directory, transverse_cells=()):
    base = Path(base_directory)
    shards = []
    for path in shard_paths:
        path = Path(path)
        relative = path.relative_to(base)
        shards.append(LedgerShard(relative.as_posix(), file_sha256(path)))
    return CampaignManifest(tuple(tuple(axis) for axis in domain), symmetry,
                            int(symmetry_factor), required_grade, tuple(shards),
                            tuple(transverse_cells))
