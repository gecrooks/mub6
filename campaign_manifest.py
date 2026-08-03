"""Canonical campaign manifest suitable for an independent rechecker."""

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

from certificate_result import CertificateGrade
from ledger_bits import bits_float, float_bits


SCHEMA = "mub6-campaign-manifest-v1"


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
class CampaignManifest:
    domain: tuple[tuple[float, float], ...]
    symmetry: str
    symmetry_factor: int
    required_grade: CertificateGrade
    ledger_shards: tuple[LedgerShard, ...]

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

    def _unsigned_dict(self):
        return {
            "schema": SCHEMA,
            "domain_bits": [[float_bits(lo), float_bits(hi)]
                            for lo, hi in self.domain],
            "symmetry": self.symmetry,
            "symmetry_factor": self.symmetry_factor,
            "required_grade": self.required_grade.name,
            "ledger_shards": [shard.as_dict() for shard in self.ledger_shards],
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
                   shard_paths, *, base_directory):
    base = Path(base_directory)
    shards = []
    for path in shard_paths:
        path = Path(path)
        relative = path.relative_to(base)
        shards.append(LedgerShard(relative.as_posix(), file_sha256(path)))
    return CampaignManifest(tuple(tuple(axis) for axis in domain), symmetry,
                            int(symmetry_factor), required_grade, tuple(shards))
