"""Download and verify the ChEMBL SQLite distribution.

The data-source version is part of the reproducibility contract, exactly like
the RDKit version: ChEMBL releases differ in their activity records and in the
coverage of `predicted_binding_domains`, so a result is only reproducible if
the release is stated. This script therefore refuses to download "latest"
without recording which numbered release that resolved to, and it writes the
resolved version plus the checksum into docs/versions.md.

Usage
-----
List the releases available on the EBI FTP mirror::

    python scripts/00_download_chembl.py --list

Download a specific release (recommended -- pin it, do not track "latest")::

    python scripts/00_download_chembl.py --version 35

Verify an archive that is already on disk and extract it::

    python scripts/00_download_chembl.py --version 35 --skip-download

Disk budget: the archive is roughly 6 GB and expands to roughly 25 GB. Both
are written under data/raw/chembl/.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import sys
import tarfile
import threading
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path

FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb"
REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = REPO_ROOT / "data" / "raw" / "chembl"
VERSIONS_DOC = REPO_ROOT / "docs" / "versions.md"

CHUNK = 1 << 20  # 1 MiB


def fetch_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode("utf-8", errors="replace")


def list_releases() -> list[str]:
    """Return the numbered releases visible in the FTP directory listing."""
    listing = fetch_text(f"{FTP_BASE}/releases/")
    versions = sorted({m for m in re.findall(r"chembl_(\d+)/", listing)}, key=int)
    return versions


class Progress:
    """Thread-safe byte counter that rewrites a single progress line."""

    def __init__(self, total: int, already_have: int = 0) -> None:
        self.total = total
        self.done = already_have
        self._lock = threading.Lock()

    def add(self, n: int) -> None:
        with self._lock:
            self.done += n
            if self.total:
                pct = 100 * self.done / self.total
                print(
                    f"\r  {self.done / 1e9:.2f} / {self.total / 1e9:.2f} GB "
                    f"({pct:.1f}%)",
                    end="",
                    flush=True,
                )


def probe(url: str) -> tuple[int, bool]:
    """Return (content length, whether the server honours byte ranges)."""
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=60) as response:
        total = int(response.headers.get("Content-Length", 0))
        accepts = response.headers.get("Accept-Ranges", "none").lower() == "bytes"
    return total, accepts


def fetch_range(
    url: str, start: int, end: int, out_path: Path, progress: Progress
) -> None:
    """Download one byte range [start, end] inclusive into its own file."""
    request = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
    with urllib.request.urlopen(request, timeout=300) as response:
        with out_path.open("wb") as handle:
            while True:
                chunk = response.read(CHUNK)
                if not chunk:
                    break
                handle.write(chunk)
                progress.add(len(chunk))

    expected = end - start + 1
    actual = out_path.stat().st_size
    if actual != expected:
        raise OSError(f"range {start}-{end}: got {actual} bytes, expected {expected}")


def download(url: str, destination: Path, connections: int = 8) -> None:
    """Download `url` to `destination`, resuming and parallelising if possible.

    EBI throttles per connection -- a single stream caps out around 1.5 MB/s
    regardless of local bandwidth -- but it serves byte ranges, so the remaining
    bytes are split across several connections. An interrupted download is
    resumed from the `.part` file rather than restarted: only whole, verified
    chunks are ever appended to it.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")

    total, accepts_ranges = probe(url)
    have = partial.stat().st_size if partial.exists() else 0

    if total and have >= total:
        print(f"  already complete ({have / 1e9:.2f} GB)")
        partial.replace(destination)
        return

    if have:
        print(f"  resuming from {have / 1e9:.2f} GB")

    # Fall back to one stream when ranges are unsupported or the size is unknown.
    if not accepts_ranges or not total or connections < 2:
        if not accepts_ranges and have:
            print("  server does not honour ranges; restarting from zero")
            partial.unlink()
            have = 0
        progress = Progress(total, have)
        headers = {"Range": f"bytes={have}-"} if have else {}
        request = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(request, timeout=300) as response:
            with partial.open("ab" if have else "wb") as handle:
                while True:
                    chunk = response.read(CHUNK)
                    if not chunk:
                        break
                    handle.write(chunk)
                    progress.add(len(chunk))
        print()
        partial.replace(destination)
        return

    # Split the remaining bytes into one contiguous range per connection.
    remaining = total - have
    span = remaining // connections
    ranges = []
    for i in range(connections):
        start = have + i * span
        end = total - 1 if i == connections - 1 else start + span - 1
        ranges.append((i, start, end))

    print(f"  {remaining / 1e9:.2f} GB left over {connections} connections")
    progress = Progress(total, have)
    chunk_dir = partial.parent / f"{partial.name}.chunks"
    chunk_dir.mkdir(exist_ok=True)
    chunk_paths = {i: chunk_dir / f"{i:03d}" for i, _, _ in ranges}

    with ThreadPoolExecutor(max_workers=connections) as pool:
        futures = {
            pool.submit(fetch_range, url, start, end, chunk_paths[i], progress): i
            for i, start, end in ranges
        }
        for future in as_completed(futures):
            # Re-raise the first failure; the .part file is untouched, so the
            # next run resumes from where this one started.
            future.result()
    print()

    # All ranges succeeded: append them in order.
    with partial.open("ab") as handle:
        for i, _, _ in ranges:
            with chunk_paths[i].open("rb") as chunk_file:
                shutil.copyfileobj(chunk_file, handle, length=CHUNK)
    shutil.rmtree(chunk_dir)

    final = partial.stat().st_size
    if final != total:
        raise OSError(f"assembled {final} bytes, expected {total}")
    partial.replace(destination)


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expected_checksum(version: str, archive_name: str) -> str | None:
    """Read the published checksum for this archive, if the mirror offers one."""
    for name in ("checksums.txt", "CHECKSUMS"):
        url = f"{FTP_BASE}/releases/chembl_{version}/{name}"
        try:
            text = fetch_text(url)
        except Exception:
            continue
        for line in text.splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[-1].lstrip("*").endswith(archive_name):
                return parts[0].lower()
    return None


def extract(archive: Path, target_db: Path) -> None:
    """Extract the single .db member of the archive to `target_db`.

    The archive is an untrusted download, so members are not extracted
    wholesale: we locate the database file and copy only that stream out.
    """
    target_db.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as tar:
        members = [m for m in tar.getmembers() if m.isfile() and m.name.endswith(".db")]
        if len(members) != 1:
            raise SystemExit(
                f"expected exactly one .db member, found {len(members)}: "
                f"{[m.name for m in members]}"
            )
        member = members[0]
        print(f"  extracting {member.name} ({member.size / 1e9:.2f} GB)")
        source = tar.extractfile(member)
        if source is None:
            raise SystemExit(f"could not read {member.name} from the archive")
        with source, target_db.open("wb") as handle:
            shutil.copyfileobj(source, handle, length=CHUNK)


def record_version(version: str, archive_name: str, checksum: str) -> None:
    """Append the resolved release and checksum to docs/versions.md."""
    VERSIONS_DOC.parent.mkdir(parents=True, exist_ok=True)
    entry = (
        f"\n## ChEMBL\n\n"
        f"- Release / Sürüm: **ChEMBL {version}**\n"
        f"- Archive / Arşiv: `{archive_name}`\n"
        f"- SHA-256: `{checksum}`\n"
        f"- Downloaded / İndirilme tarihi: {date.today().isoformat()}\n"
        f"- Source / Kaynak: {FTP_BASE}/releases/chembl_{version}/\n"
    )
    with VERSIONS_DOC.open("a", encoding="utf-8") as handle:
        handle.write(entry)
    print(f"  recorded in {VERSIONS_DOC.relative_to(REPO_ROOT)}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="list available releases")
    parser.add_argument("--version", help="release number to download, e.g. 35")
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="verify and extract an archive already present under data/raw/chembl/",
    )
    parser.add_argument(
        "--keep-archive",
        action="store_true",
        help="keep the .tar.gz after extraction (needs ~6 GB extra)",
    )
    parser.add_argument(
        "--connections",
        type=int,
        default=8,
        help="parallel connections (default 8; EBI throttles each stream to "
             "about 1.5 MB/s, so one connection is the slow path). Use 1 to "
             "force a single stream.",
    )
    args = parser.parse_args()

    if args.list:
        versions = list_releases()
        print("Available ChEMBL releases:")
        print("  " + ", ".join(versions))
        print(f"\nMost recent: {versions[-1]}")
        print("Pin one explicitly:  --version <number>")
        return 0

    if not args.version:
        parser.error("either --list or --version is required")

    version = args.version.lstrip("chembl_")
    archive_name = f"chembl_{version}_sqlite.tar.gz"
    archive = RAW_DIR / archive_name
    target_db = RAW_DIR / "chembl.db"
    url = f"{FTP_BASE}/releases/chembl_{version}/{archive_name}"

    if target_db.exists():
        print(f"{target_db} already exists; nothing to do.")
        print("Delete it first if you intend to switch releases.")
        return 0

    if not args.skip_download:
        print(f"Downloading {url}")
        print("  (5.8 GB; interrupting is safe, the next run resumes)")
        download(url, archive, connections=args.connections)
    elif not archive.exists():
        raise SystemExit(f"--skip-download given but {archive} is missing")

    print("Computing SHA-256")
    checksum = sha256sum(archive)
    print(f"  {checksum}")

    published = expected_checksum(version, archive_name)
    if published is None:
        print("  WARNING: no published checksum found on the mirror.")
        print("  The archive is recorded as-is; verify manually if it matters.")
    elif published != checksum:
        raise SystemExit(
            f"checksum mismatch\n  published: {published}\n  local:     {checksum}\n"
            "The download is corrupt or the release was republished. Delete the "
            "archive and retry."
        )
    else:
        print("  matches the published checksum")

    print("Extracting")
    extract(archive, target_db)

    if not args.keep_archive:
        archive.unlink()
        print(f"  removed {archive.name}")

    record_version(version, archive_name, checksum)
    size_gb = target_db.stat().st_size / 1e9
    print(f"\nDone. {target_db.relative_to(REPO_ROOT)} ({size_gb:.1f} GB)")
    print("Next: python scripts/01_inspect_chembl.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
