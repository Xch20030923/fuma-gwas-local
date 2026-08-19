#!/usr/bin/env python3
"""Fixed-path FUMA-compatible GWAS detector, queue, runner, and interface.

The script deliberately uses only the Python standard library for orchestration.
The actual SNP2GENE production work is delegated to the audited local project
launcher, so reference handling and result semantics stay in one place.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from typing import Any, Iterable, Iterator, TextIO
from xml.sax.saxutils import escape as xml_escape


PROJECT_ROOT = Path("/media/desk16/iy19619/FUMA-compatible-SNP2GENE")
POSTGWAS_ROOT = Path("/media/desk16/iy19619/iyun8003/post-GWAS分析")
POSTGWAS_CODE_ROOT = POSTGWAS_ROOT / "post-GWAS代码文件"
LOCAL_ROOT = POSTGWAS_ROOT / "FUMA_local_runs"
QUEUE_ROOT = LOCAL_ROOT / "queue"
JOBS_ROOT = LOCAL_ROOT / "jobs"
PROFILE = PROJECT_ROOT / "config/profiles/fuma_v2.1.6_eur_hg19_ensembl_v102.yaml"
PRODUCTION_LAUNCHER = PROJECT_ROOT / "scripts/build_candidate_annovar_background.sh"

TEXT_SUFFIXES = {".txt", ".tsv", ".csv", ".tab", ".gz", ".zip"}
EXCLUDED_DIR_NAMES = {
    ".git",
    ".codex",
    "FUMA_local_runs",
    "FUMA_ready",
    "FUMA_zips",
    "runs",
    "results",
    "output",
    "outputs",
}
FUMA_OUTPUT_NAMES = {
    "snps.txt",
    "snps_positional.txt",
    "annot.txt",
    "annov.txt",
    "genes.txt",
    "genes_positional.txt",
    "GenomicRiskLoci.txt",
    "IndSigSNPs.txt",
    "leadSNPs.txt",
    "ld.txt",
}

ALIASES = {
    "chr": {
        "chr", "chrom", "chromosome", "chromosomeid", "chromosomecode",
        "chrnum", "chromnum", "chromosome_number",
    },
    "pos": {
        "pos", "bp", "position", "basepair", "basepairlocation",
        "basepairposition", "physicalposition", "location", "start",
    },
    "snp": {
        "snp", "rsid", "rsnumber", "marker", "markername", "variant",
        "variantid", "variantidentifier", "id", "snpid", "snpname",
    },
    "p": {
        "p", "pvalue", "pval", "pvalueoverall", "pplaco", "placop",
        "placopvalue", "metap", "pmeta", "pwald", "p_gwas", "p_assoc",
    },
    "neglogp": {
        "neglog10p", "minuslog10p", "mlogp", "log10p", "neglogp",
        "minuslogp", "logp",
    },
}


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def json_dump(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.strip().lower())


def clean_chr(value: str) -> str | None:
    value = value.strip().replace("\"", "")
    if value.lower().startswith("chr"):
        value = value[3:]
    if value in {"23", "x", "X"}:
        return "X"
    if value.isdigit() and 1 <= int(value) <= 22:
        return str(int(value))
    return None


def parse_coord(value: str) -> tuple[str, int] | None:
    match = re.match(r"^(?:chr)?([0-9]{1,2}|X|x)[:_\-]([0-9]+)", value.strip())
    if not match:
        return None
    chrom = clean_chr(match.group(1))
    try:
        pos = int(match.group(2))
    except ValueError:
        return None
    return (chrom, pos) if chrom and pos > 0 else None


def is_missing(value: str) -> bool:
    return not value.strip() or value.strip().upper() in {"NA", "N/A", ".", "NULL", "NAN"}


def split_fields(line: str, delimiter: str) -> list[str]:
    if delimiter == "whitespace":
        return line.strip().split()
    return next(csv.reader([line], delimiter=delimiter))


def choose_delimiter(line: str) -> str:
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    if ";" in line and line.count(";") >= 2:
        return ";"
    return "whitespace"


def iter_lines(path: Path, member: str | None = None) -> Iterator[str]:
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            chosen = member
            if chosen is None:
                names = [n for n in archive.namelist() if not n.endswith("/")]
                chosen = next((n for n in names if not Path(n).name.startswith(".")), names[0])
            with archive.open(chosen, "r") as raw:
                with TextIOWrapper(raw, encoding="utf-8-sig", errors="replace", newline="") as handle:
                    yield from handle
        return
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as handle:
        yield from handle


def text_members(path: Path) -> list[str | None]:
    if path.suffix.lower() != ".zip":
        return [None]
    with zipfile.ZipFile(path) as archive:
        return [
            name for name in archive.namelist()
            if not name.endswith("/") and Path(name).suffix.lower() in TEXT_SUFFIXES | {".bgz"}
        ]


def sample_content(path: Path, member: str | None, limit: int = 60) -> tuple[list[tuple[int, str]], str]:
    lines: list[tuple[int, str]] = []
    for number, line in enumerate(iter_lines(path, member), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        lines.append((number, stripped))
        if len(lines) >= limit:
            break
    return lines, choose_delimiter(lines[0][1]) if lines else "whitespace"


def alias_kind(name: str) -> str | None:
    normalized = normalize_name(name)
    for kind, values in ALIASES.items():
        if normalized in {normalize_name(v) for v in values}:
            return kind
    return None


def looks_like_header(fields: list[str]) -> bool:
    kinds = {alias_kind(field) for field in fields}
    kinds.discard(None)
    return len(kinds & {"chr", "pos", "snp", "p", "neglogp"}) >= 2


def numeric(value: str) -> float | None:
    try:
        parsed = float(value.strip())
    except (AttributeError, TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def infer_indexes(rows: list[list[str]]) -> dict[str, int | None]:
    if not rows:
        return {"chr": None, "pos": None, "snp": None, "p": None, "neglogp": None}
    width = max(len(row) for row in rows)
    columns = [[row[i] if i < len(row) else "" for row in rows] for i in range(width)]
    snp_scores = [sum(bool(re.match(r"^(?:rs[0-9]+|(?:chr)?[0-9XY]+[:_\-][0-9]+)", x, re.I)) for x in col) for col in columns]
    chr_scores = [sum(clean_chr(x) is not None for x in col) for col in columns]
    pos_scores = [sum(bool(x.isdigit()) and int(x) > 0 for x in col if x.strip()) for col in columns]
    p_scores = [sum((numeric(x) is not None and 0 <= float(x) <= 1) for x in col) for col in columns]
    snp_index = max(range(width), key=lambda i: snp_scores[i]) if width else None
    chr_index = max(range(width), key=lambda i: chr_scores[i]) if width else None
    pos_index = max(range(width), key=lambda i: pos_scores[i]) if width else None
    p_index = max(range(width), key=lambda i: p_scores[i]) if width else None
    if snp_index is not None and snp_scores[snp_index] == 0:
        snp_index = None
    if chr_index is not None and chr_scores[chr_index] == 0:
        chr_index = None
    if pos_index is not None and pos_scores[pos_index] == 0:
        pos_index = None
    if p_index is not None and p_scores[p_index] == 0:
        p_index = None
    return {"chr": chr_index, "pos": pos_index, "snp": snp_index, "p": p_index, "neglogp": None}


def detect_member(path: Path, member: str | None) -> dict[str, Any]:
    lines, delimiter = sample_content(path, member)
    if not lines:
        return {"status": "FAIL_EMPTY", "path": str(path), "member": member}
    first_fields = split_fields(lines[0][1], delimiter)
    header = looks_like_header(first_fields)
    if header:
        indexes: dict[str, int | None] = {"chr": None, "pos": None, "snp": None, "p": None, "neglogp": None}
        for idx, field in enumerate(first_fields):
            kind = alias_kind(field)
            if kind and indexes[kind] is None:
                indexes[kind] = idx
        data_rows = [split_fields(line, delimiter) for _, line in lines[1:]]
        header_line = lines[0][0]
        header_fields = first_fields
    else:
        data_rows = [split_fields(line, delimiter) for _, line in lines]
        indexes = infer_indexes(data_rows)
        header_line = None
        header_fields = []
    required = {"chr", "pos", "p"} <= {k for k, v in indexes.items() if v is not None} and (
        indexes["snp"] is not None or True
    )
    score = sum(indexes[k] is not None for k in ("chr", "pos", "snp", "p", "neglogp"))
    status = "PASS" if score >= 3 and (indexes["p"] is not None or indexes["neglogp"] is not None) and (indexes["chr"] is not None or indexes["snp"] is not None) else "REVIEW"
    return {
        "status": status,
        "path": str(path),
        "member": member,
        "delimiter": delimiter,
        "header_detected": header,
        "header_line": header_line,
        "header": header_fields,
        "indexes": indexes,
        "score": score,
        "sample_rows": len(data_rows),
        "sample_columns": max((len(row) for row in data_rows), default=0),
        "required_shape_detected": required,
    }


def detect_source(path: Path) -> dict[str, Any]:
    members = text_members(path)
    detections = [detect_member(path, member) for member in members]
    detections.sort(key=lambda item: (item.get("status") == "PASS", item.get("score", 0), item.get("sample_rows", 0)), reverse=True)
    result = detections[0] if detections else {"status": "FAIL_NO_TEXT_MEMBER"}
    result["source_sha256"] = sha256(path)
    result["source_size_bytes"] = path.stat().st_size
    result["source_path"] = str(path.resolve())
    if len(detections) > 1:
        result["members_considered"] = detections
    return result


def candidate_files(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    if not path.is_dir():
        return []
    found: list[Path] = []
    for current, dirs, files in os.walk(path):
        dirs[:] = [d for d in dirs if d not in EXCLUDED_DIR_NAMES and not d.startswith(".")]
        for filename in files:
            candidate = Path(current) / filename
            if candidate.suffix.lower() not in TEXT_SUFFIXES:
                continue
            if candidate.name in FUMA_OUTPUT_NAMES or candidate.name.endswith(".status.json"):
                continue
            found.append(candidate)
    return sorted(found)


def expand_inputs(paths: Iterable[str]) -> list[Path]:
    expanded: list[Path] = []
    for value in paths:
        path = Path(value).expanduser().resolve()
        if path.is_file():
            expanded.append(path)
        elif path.is_dir():
            expanded.extend(candidate_files(path))
        else:
            raise FileNotFoundError(path)
    unique: list[Path] = []
    seen: set[Path] = set()
    for path in expanded:
        if path not in seen:
            unique.append(path)
            seen.add(path)
    return unique


def value_at(row: list[str], index: int | None) -> str:
    return row[index].strip() if index is not None and index < len(row) else ""


def normalize_input(source: Path, output: Path, detection: dict[str, Any]) -> dict[str, Any]:
    delimiter = detection.get("delimiter", "whitespace")
    indexes = detection.get("indexes", {})
    member = detection.get("member")
    header_line = detection.get("header_line")
    rows_seen = 0
    rows_kept = 0
    dropped: dict[str, int] = {}
    records: dict[tuple[str, int, str], tuple[float, str]] = {}

    def drop(reason: str) -> None:
        dropped[reason] = dropped.get(reason, 0) + 1

    for number, line in enumerate(iter_lines(source, member), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or (header_line and number == header_line):
            continue
        fields = split_fields(stripped, delimiter)
        rows_seen += 1
        snp = value_at(fields, indexes.get("snp"))
        chrom = clean_chr(value_at(fields, indexes.get("chr")))
        position_value = value_at(fields, indexes.get("pos"))
        if (chrom is None or not position_value) and snp:
            coord = parse_coord(snp)
            if coord:
                chrom, position_value = coord[0], str(coord[1])
        if chrom is None:
            drop("invalid_chr")
            continue
        try:
            position = int(float(position_value))
        except (TypeError, ValueError):
            drop("invalid_pos")
            continue
        if position <= 0:
            drop("invalid_pos")
            continue
        p_raw = value_at(fields, indexes.get("p"))
        transform = "p"
        if indexes.get("neglogp") is not None:
            p_raw = value_at(fields, indexes.get("neglogp"))
            transform = "neglog10p"
        p_value = numeric(p_raw)
        if p_value is None:
            drop("invalid_p")
            continue
        if transform == "neglog10p":
            if p_value < 0:
                drop("invalid_neglog10p")
                continue
            p_value = 10.0 ** (-min(p_value, 300.0))
        elif p_value == 0:
            p_value = 1e-300
        if not (0 < p_value <= 1):
            drop("invalid_p")
            continue
        if is_missing(snp):
            snp = f"{chrom}:{position}"
        key = (chrom, position, snp)
        previous = records.get(key)
        if previous is None or p_value < previous[0]:
            records[key] = (p_value, snp)
        rows_kept += 1

    duplicates = rows_kept - len(records)
    if duplicates:
        dropped["duplicate_variant_best_p_retained"] = duplicates
    output.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(records.items(), key=lambda item: (item[0][0] != "X", int(item[0][0]) if item[0][0] != "X" else 23, item[0][1], item[0][2]))
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chr", "pos", "SNP", "p.placo"])
        for (chrom, position, _), (p_value, snp) in ordered:
            writer.writerow([chrom, position, snp, f"{p_value:.12g}"])
    return {
        "status": "PASS" if records else "FAIL_NO_VALID_ROWS",
        "source_path": str(source.resolve()),
        "source_member": member,
        "output_path": str(output.resolve()),
        "rows_seen": rows_seen,
        "rows_valid_before_dedup": rows_kept,
        "rows_written": len(records),
        "rows_dropped": rows_seen - rows_kept,
        "drop_reasons": dict(sorted(dropped.items())),
        "columns": {key: indexes.get(key) for key in ("chr", "pos", "snp", "p", "neglogp")},
        "output_sha256": sha256(output) if output.exists() else None,
    }


def safe_id(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    return value[:90]


def unique_run_id(source: Path, source_hash: str, requested: str | None = None) -> str:
    if requested:
        base = safe_id(requested)
    else:
        stem = source.name
        for suffix in (".gz", ".bgz", ".zip"):
            if stem.lower().endswith(suffix):
                stem = stem[: -len(suffix)]
        base = safe_id(Path(stem).stem) or "gwas"
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    candidate = f"{base}_{stamp}_{source_hash[:8]}"
    if requested:
        candidate = safe_id(requested)
    result_path = PROJECT_ROOT / "runs" / candidate
    if not result_path.exists() and not (JOBS_ROOT / f"{candidate}.json").exists():
        return candidate
    index = 2
    while (PROJECT_ROOT / "runs" / f"{candidate}_{index}").exists() or (JOBS_ROOT / f"{candidate}_{index}.json").exists():
        index += 1
    return f"{candidate}_{index}"


def job_path(run_id: str) -> Path:
    return JOBS_ROOT / f"{run_id}.json"


def read_job(run_id: str) -> dict[str, Any]:
    path = job_path(run_id)
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def write_job(job: dict[str, Any]) -> None:
    job["updated_at_utc"] = now_iso()
    json_dump(job_path(job["run_id"]), job)


def process_exists(pid: int | None) -> bool:
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def scheduler_alive(queue_root: Path) -> bool:
    pid_file = queue_root / "scheduler.pid"
    try:
        pid = int(pid_file.read_text(encoding="utf-8").strip())
    except (OSError, ValueError):
        return False
    return process_exists(pid)


def start_scheduler(workers: int, queue_root: Path) -> dict[str, Any]:
    queue_root.mkdir(parents=True, exist_ok=True)
    log_path = queue_root / "scheduler.log"
    if scheduler_alive(queue_root):
        return {"status": "ALREADY_RUNNING", "pid_file": str(queue_root / "scheduler.pid"), "log": str(log_path)}
    command = [sys.executable, str(Path(__file__).resolve()), "worker", "--workers", str(workers), "--queue-root", str(queue_root)]
    with log_path.open("a", encoding="utf-8") as log_handle:
        process = subprocess.Popen(
            command,
            cwd=str(PROJECT_ROOT),
            stdin=subprocess.DEVNULL,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    (queue_root / "scheduler.pid").write_text(str(process.pid) + "\n", encoding="utf-8")
    return {"status": "STARTED", "pid": process.pid, "pid_file": str(queue_root / "scheduler.pid"), "log": str(log_path)}


def discover_jobs(states: set[str]) -> list[dict[str, Any]]:
    jobs: list[dict[str, Any]] = []
    for path in sorted(JOBS_ROOT.glob("*.json")):
        try:
            job = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if job.get("state") in states:
            jobs.append(job)
    return jobs


def find_raw_results(run_id: str) -> Path:
    path = PROJECT_ROOT / "runs" / run_id / "results"
    if not path.is_dir():
        raise FileNotFoundError(path)
    return path


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        first = handle.readline()
        if not first:
            return []
        delimiter = "\t" if "\t" in first else ","
        handle.seek(0)
        return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def write_table(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fieldnames} for row in rows)


def unique_values(values: Iterable[str]) -> list[str]:
    output: list[str] = []
    seen: set[str] = set()
    for raw in values:
        value = str(raw or "").strip()
        if not value or value.upper() in {"NA", "N/A", "NAN", "NULL"}:
            continue
        if value not in seen:
            output.append(value)
            seen.add(value)
    return output


def symlink_result_files(raw: Path, interface_dir: Path) -> dict[str, Any]:
    linked: list[str] = []
    missing: list[str] = []
    for name in sorted(FUMA_OUTPUT_NAMES | {"annov.stats.txt", "EUR.annov.count", "gwascatalog.txt", "gwascatalog.status.json", "params.config", "README", "reference_manifest.tsv", "input_qc.tsv"}):
        source = raw / name
        destination = interface_dir / name
        if source.is_file():
            if destination.exists() or destination.is_symlink():
                destination.unlink()
            destination.symlink_to(os.path.relpath(source, interface_dir))
            linked.append(name)
        elif name in {"snps.txt", "genes.txt", "annov.txt", "GenomicRiskLoci.txt"}:
            missing.append(name)
    return {"linked": linked, "missing_required": missing}


def write_xlsx(path: Path, sheets: list[tuple[str, list[dict[str, Any]]]]) -> None:
    """Write a dependency-free XLSX with inline strings."""
    def col_name(index: int) -> str:
        result = ""
        index += 1
        while index:
            index, remainder = divmod(index - 1, 26)
            result = chr(65 + remainder) + result
        return result

    def sheet_xml(rows: list[dict[str, Any]]) -> str:
        fields = list(rows[0].keys()) if rows else []
        all_rows = [fields] + [[row.get(field, "") for field in fields] for row in rows]
        body: list[str] = []
        for row_index, values in enumerate(all_rows, start=1):
            cells: list[str] = []
            for col_index, value in enumerate(values):
                text = xml_escape(str(value if value is not None else ""))
                cells.append(f'<c r="{col_name(col_index)}{row_index}" t="inlineStr"><is><t>{text}</t></is></c>')
            body.append(f'<row r="{row_index}">' + "".join(cells) + "</row>")
        return '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>' + "".join(body) + "</sheetData></worksheet>"

    path.parent.mkdir(parents=True, exist_ok=True)
    content_types = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">', '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>', '<Default Extension="xml" ContentType="application/xml"/>', '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>']
    rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>']
    workbook_sheets: list[str] = []
    workbook_rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for index, (name, rows) in enumerate(sheets, start=1):
            archive.writestr(f"xl/worksheets/sheet{index}.xml", sheet_xml(rows))
            content_types.append(f'<Override PartName="/xl/worksheets/sheet{index}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
            workbook_sheets.append(f'<sheet name="{xml_escape(name[:31])}" sheetId="{index}" r:id="rId{index}"/>')
            workbook_rels.append(f'<Relationship Id="rId{index}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{index}.xml"/>')
        workbook_rels.append("</Relationships>")
        archive.writestr("[Content_Types].xml", "".join(content_types) + "</Types>")
        archive.writestr("_rels/.rels", "".join(rels))
        archive.writestr("xl/_rels/workbook.xml.rels", "".join(workbook_rels))
        archive.writestr("xl/workbook.xml", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>' + "".join(workbook_sheets) + "</sheets></workbook>")


def interface_run(job: dict[str, Any], coloc_path: Path | None = None) -> dict[str, Any]:
    raw = find_raw_results(job["run_id"])
    interface_dir = Path(job["job_root"]) / "interface"
    interface_dir.mkdir(parents=True, exist_ok=True)
    links = symlink_result_files(raw, interface_dir)
    loci = read_table(raw / "GenomicRiskLoci.txt")
    snps = read_table(raw / "snps.txt")
    genes = read_table(raw / "genes.txt")
    annov = read_table(raw / "annov.txt")
    nearest_by_snp: dict[str, list[str]] = {}
    for row in snps:
        nearest_by_snp.setdefault(row.get("rsID", ""), []).append(row.get("nearestGene", ""))
    locus_rows: list[dict[str, Any]] = []
    for locus in loci:
        lead_genes: list[str] = []
        for rsid in str(locus.get("LeadSNPs", "")).split(";"):
            lead_genes.extend(nearest_by_snp.get(rsid, []))
        result = dict(locus)
        chrom = result.get("chr", "")
        result["CHR"] = f"chr{chrom}" if chrom else ""
        result["cytoBand"] = "NA"
        result["cytoBand_source"] = "not bundled; use QTLMR data_info_cytoBand(build=37) when required"
        result["LeadSNP_NearestGene"] = ";".join(unique_values(lead_genes)) or "NA"
        locus_rows.append(result)
    coloc_rows = read_table(coloc_path) if coloc_path and coloc_path.is_file() else []
    if coloc_rows:
        by_risk = {str(row.get("risk_locus", row.get("GenomicLocus", ""))): row for row in coloc_rows}
        coloc_fields = unique_values(key for row in coloc_rows for key in row.keys())
        for result in locus_rows:
            match = by_risk.get(str(result.get("GenomicLocus", "")), {})
            for key in coloc_fields:
                if key in {"risk_locus", "GenomicLocus"}:
                    continue
                result[f"coloc_{key}"] = match.get(key, "")
            result["coloc_status"] = "matched" if match else "not_matched"
    else:
        for result in locus_rows:
            result["coloc_status"] = "not_supplied"
    locus_csv = interface_dir / "GenomicRiskLoci_with_NearestGene_and_Coloc.csv"
    write_table(locus_csv, locus_rows)
    write_xlsx(interface_dir / "GenomicRiskLoci_with_NearestGene_and_Coloc.xlsx", [("GenomicRiskLoci+coloc", locus_rows)])

    nearest = unique_values(
        gene
        for row in locus_rows
        for value in str(row.get("LeadSNP_NearestGene", "")).split(";")
        for gene in re.split(r"[:;]", value)
    )
    v2g = unique_values(row.get("symbol", "") for row in genes)
    exonic = unique_values(row.get("symbol", "") for row in annov if str(row.get("annot", "")).lower() == "exonic")
    size = max(len(nearest), len(v2g), len(exonic), 1)
    tri_rows = [{"Nearest": nearest[i] if i < len(nearest) else "", "V2G": v2g[i] if i < len(v2g) else "", "Exonic": exonic[i] if i < len(exonic) else ""} for i in range(size)]
    write_table(interface_dir / "FUMA_V2G_Nearest_Exonic.csv", tri_rows, ["Nearest", "V2G", "Exonic"])
    final_xlsx = interface_dir / f"{job['run_id']}_Final.xlsx"
    write_xlsx(final_xlsx, [("Nearest_V2G_Exonic", tri_rows), ("GenomicRiskLoci+coloc", locus_rows)])

    manifest_files: list[dict[str, Any]] = []
    for path in sorted(interface_dir.iterdir()):
        if path.is_file() and not path.is_symlink():
            manifest_files.append({"name": path.name, "path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    manifest = {
        "status": "PASS" if not links["missing_required"] else "FAIL_MISSING_REQUIRED_FUMA_FILES",
        "run_id": job["run_id"],
        "raw_results": str(raw),
        "interface_dir": str(interface_dir),
        "postgwas_code_root": str(POSTGWAS_CODE_ROOT),
        "linked_raw_files": links,
        "derived_files": manifest_files,
        "coloc_file": str(coloc_path) if coloc_path else None,
        "downstream_contract": {
            "nearest_coloc_csv": str(locus_csv),
            "nearest_coloc_xlsx": str(interface_dir / "GenomicRiskLoci_with_NearestGene_and_Coloc.xlsx"),
            "final_xlsx": str(final_xlsx),
            "raw_fuma_directory_can_be_passed_to_existing_R_function": True,
        },
        "cytoband_note": "The compatible CSV includes cytoBand=NA; the existing R workflow can fill it with QTLMR::data_info_cytoBand(build=37).",
        "created_at_utc": now_iso(),
    }
    json_dump(interface_dir / "FUMA_postGWAS_manifest.json", manifest)
    return manifest


def build_report(job: dict[str, Any], interface: dict[str, Any] | None = None) -> dict[str, Any]:
    raw = PROJECT_ROOT / "runs" / job["run_id"] / "results"
    required = ["IndSigSNPs.txt", "leadSNPs.txt", "GenomicRiskLoci.txt", "snps.txt", "ld.txt", "genes.txt", "annov.txt", "annot.txt"]
    present = [name for name in required if (raw / name).is_file()]
    report = {
        "status": "PASS_PRACTICAL_CANDIDATE" if len(present) == len(required) else "INCOMPLETE",
        "fidelity_level": "practical_candidate",
        "strict_1_to_1_usable": False,
        "run_id": job["run_id"],
        "input": job.get("input", {}),
        "profile": str(PROFILE),
        "official_fuma_commit": "9ccff60570ea06a43bd7fa77aeb62920ad271df4",
        "required_key_outputs": required,
        "present_key_outputs": present,
        "reference_limitations": [
            "FUMA private dbSNP146/RsMerge146 processed snapshot is unavailable; public structural candidates are used and flagged.",
            "FUMA private EUR precomputed LD/frequency archive is unavailable; the local production backend uses the audited EUR PLINK reference and dynamic LD.",
            "Processed private annotation snapshots are unavailable; annotation values may differ while core locus membership remains the priority result.",
        ],
        "runtime_estimate": {
            "input_detection_and_normalization": "seconds to minutes, mainly proportional to compressed input size",
            "warm_single_pair": "approximately 10-30 minutes for a typical PLACO-sized input after all references are ready; verify from job timestamps",
            "cold_reference_preparation": "hours to days and should be treated as a one-time preparation, not a per-GWAS runtime",
            "parallel_recommendation": "workers=2 by default; workers=3-4 only after memory and PLINK I/O are checked",
        },
        "parallel_jobs_supported": True,
        "postgwas_interface": interface,
        "created_at_utc": now_iso(),
    }
    report_path = Path(job["job_root"]) / "reproducibility_report.json"
    json_dump(report_path, report)
    markdown = [
        f"# FUMA local reproducibility report: {job['run_id']}",
        "",
        f"- Status: `{report['status']}`",
        f"- Fidelity: `{report['fidelity_level']}`; strict 1:1 usable: `{report['strict_1_to_1_usable']}`",
        f"- Input rows written: `{job.get('normalization', {}).get('rows_written', 'NA')}`",
        "",
        "## Required outputs",
        "",
        *[f"- `{name}`: {'present' if name in present else 'MISSING'}" for name in required],
        "",
        "## Reproducibility limits",
        "",
        *[f"- {item}" for item in report["reference_limitations"]],
        "",
        "## Runtime and suspension",
        "",
        "The scheduler is detached from the calling session. It writes one job log and one JSON status file per task. The caller can stop observing immediately and later use `status` or `report`; stopping the Codex session does not stop the detached scheduler.",
    ]
    (Path(job["job_root"]) / "reproducibility_report.md").write_text("\n".join(markdown) + "\n", encoding="utf-8")
    return report


def submit_one(source: Path, requested_id: str | None, coloc: Path | None) -> dict[str, Any]:
    detection = detect_source(source)
    if detection.get("status") not in {"PASS", "REVIEW"}:
        raise ValueError(f"Unable to identify GWAS columns in {source}: {detection.get('status')}")
    run_id = unique_run_id(source, detection["source_sha256"], requested_id)
    root = LOCAL_ROOT / "jobs" / run_id
    normalized_path = root / "input" / "normalized_input.tsv"
    normalization = normalize_input(source, normalized_path, detection)
    if normalization["status"] != "PASS":
        raise ValueError(f"No valid GWAS rows after normalization: {source}")
    job = {
        "schema_version": 1,
        "run_id": run_id,
        "state": "queued",
        "created_at_utc": now_iso(),
        "updated_at_utc": now_iso(),
        "job_root": str(root),
        "source": detection,
        "input": {"normalized_path": str(normalized_path), "sha256": normalization["output_sha256"]},
        "normalization": normalization,
        "profile": str(PROFILE),
        "coloc_file": str(coloc) if coloc else None,
        "resource_mode": "production_practical_candidate",
        "scheduler": {"queue_root": str(QUEUE_ROOT), "max_workers": 2},
        "retry_count": 0,
    }
    json_dump(root / "detection.json", detection)
    json_dump(root / "normalization.json", normalization)
    write_job(job)
    return job


def cmd_detect(args: argparse.Namespace) -> int:
    sources = expand_inputs(args.paths)
    results: list[dict[str, Any]] = []
    for source in sources:
        try:
            results.append(detect_source(source))
        except Exception as exc:  # keep batch detection auditable
            results.append({"status": "FAIL_EXCEPTION", "source_path": str(source), "error": repr(exc)})
    output = Path(args.output).expanduser() if args.output else None
    if output:
        json_dump(output, {"created_at_utc": now_iso(), "results": results})
    print(json.dumps(results, indent=2, ensure_ascii=False))
    return 0 if results and all(result.get("status") == "PASS" for result in results) else 2


def cmd_submit(args: argparse.Namespace) -> int:
    sources = expand_inputs(args.paths)
    if not sources:
        raise ValueError("No candidate text files found")
    jobs: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    explicit_files = {Path(value).expanduser().resolve() for value in args.paths if Path(value).expanduser().is_file()}
    for index, source in enumerate(sources):
        detection = detect_source(source)
        if detection.get("status") != "PASS":
            skipped.append({"path": str(source), "status": detection.get("status"), "score": detection.get("score", 0), "reason": "not recognized as a GWAS table"})
            if source in explicit_files:
                raise ValueError(f"Explicit input is not recognized as a GWAS table: {source}")
            continue
        requested = args.run_id if len(sources) == 1 else (f"{args.run_id}_{index + 1}" if args.run_id else None)
        jobs.append(submit_one(source, requested, Path(args.coloc).expanduser().resolve() if args.coloc else None))
    if not jobs:
        raise ValueError("No recognizable GWAS tables were found; skipped=" + json.dumps(skipped, ensure_ascii=False))
    start = start_scheduler(max(1, min(4, args.workers)), QUEUE_ROOT)
    for job in jobs:
        job["scheduler"]["max_workers"] = max(1, min(4, args.workers))
        write_job(job)
    summary = {"status": "SUBMITTED", "jobs": jobs, "skipped": skipped, "scheduler": start, "created_at_utc": now_iso()}
    json_dump(LOCAL_ROOT / "last_submit.json", summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


def launch_job(job: dict[str, Any]) -> tuple[subprocess.Popen[str], TextIO]:
    log_path = Path(job["job_root"]) / "logs" / "fuma_pipeline.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_handle = log_path.open("a", encoding="utf-8")
    env = os.environ.copy()
    env.update({
        "FUMA_PROFILE": str(PROFILE),
        "FUMA_INPUT": job["input"]["normalized_path"],
        "FUMA_RUN_ID": job["run_id"],
        "FUMA_HISTORICAL_DIR": str(PROJECT_ROOT / ".no_historical_snapshot_for_new_run"),
    })
    command = ["bash", str(PRODUCTION_LAUNCHER), job["run_id"]]
    job["command"] = command
    job["environment"] = {key: env[key] for key in ("FUMA_PROFILE", "FUMA_INPUT", "FUMA_RUN_ID", "FUMA_HISTORICAL_DIR")}
    process = subprocess.Popen(command, cwd=str(PROJECT_ROOT), env=env, stdout=log_handle, stderr=subprocess.STDOUT, text=True)
    job["pid"] = process.pid
    job["started_at_utc"] = now_iso()
    return process, log_handle


def reconcile_running_jobs() -> None:
    for job in discover_jobs({"running"}):
        if not process_exists(job.get("pid")):
            if int(job.get("retry_count", 0)) < 1:
                job["state"] = "queued"
                job["retry_count"] = int(job.get("retry_count", 0)) + 1
                job["reconciled_at_utc"] = now_iso()
            else:
                job["state"] = "failed_unknown_process"
                job["finished_at_utc"] = now_iso()
            write_job(job)


def worker_loop(workers: int, queue_root: Path) -> int:
    queue_root.mkdir(parents=True, exist_ok=True)
    active: dict[str, tuple[subprocess.Popen[str], TextIO, dict[str, Any]]] = {}
    while True:
        reconcile_running_jobs()
        for run_id, (process, log_handle, job) in list(active.items()):
            code = process.poll()
            if code is None:
                continue
            log_handle.close()
            job["return_code"] = code
            job["finished_at_utc"] = now_iso()
            if code == 0:
                try:
                    interface = interface_run(job, Path(job["coloc_file"]) if job.get("coloc_file") else None)
                    build_report(job, interface)
                    job["state"] = "completed"
                    job["interface"] = interface
                except Exception as exc:
                    job["state"] = "postprocess_failed"
                    job["postprocess_error"] = repr(exc)
            else:
                job["state"] = "failed"
            write_job(job)
            del active[run_id]
        available = max(0, workers - len(active))
        if available:
            for job in discover_jobs({"queued"})[:available]:
                try:
                    process, log_handle = launch_job(job)
                    job["state"] = "running"
                    write_job(job)
                    active[job["run_id"]] = (process, log_handle, job)
                except Exception as exc:
                    job["state"] = "failed_to_start"
                    job["error"] = repr(exc)
                    job["finished_at_utc"] = now_iso()
                    write_job(job)
        if not active and not discover_jobs({"queued", "running"}):
            break
        time.sleep(2)
    try:
        (queue_root / "scheduler.pid").unlink()
    except FileNotFoundError:
        pass
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    jobs = [read_job(args.run_id)] if args.run_id else discover_jobs({"queued", "running", "completed", "failed", "failed_to_start", "failed_unknown_process", "postprocess_failed"})
    rows = []
    for job in jobs:
        rows.append({"run_id": job.get("run_id"), "state": job.get("state"), "rows": job.get("normalization", {}).get("rows_written"), "pid": job.get("pid", ""), "created": job.get("created_at_utc"), "finished": job.get("finished_at_utc", ""), "job_root": job.get("job_root")})
    print(json.dumps({"scheduler_alive": scheduler_alive(QUEUE_ROOT), "jobs": rows}, indent=2, ensure_ascii=False))
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    job = read_job(args.run_id)
    report_path = Path(job["job_root"]) / "reproducibility_report.json"
    if not report_path.is_file() and job.get("state") in {"completed", "postprocess_failed"}:
        report = build_report(job, job.get("interface"))
    else:
        report = json.loads(report_path.read_text(encoding="utf-8")) if report_path.is_file() else {"state": job.get("state"), "message": "Report is not ready."}
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if report.get("status", "").startswith("PASS") else 2


def cmd_interface(args: argparse.Namespace) -> int:
    job = read_job(args.run_id)
    interface = interface_run(job, Path(args.coloc).expanduser().resolve() if args.coloc else (Path(job["coloc_file"]) if job.get("coloc_file") else None))
    build_report(job, interface)
    job["interface"] = interface
    write_job(job)
    print(json.dumps(interface, indent=2, ensure_ascii=False))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Fixed-path FUMA-compatible GWAS local runner")
    sub = parser.add_subparsers(dest="command", required=True)
    detect = sub.add_parser("detect", help="detect GWAS columns without running FUMA")
    detect.add_argument("paths", nargs="+", help="GWAS files or directories")
    detect.add_argument("--output")
    detect.set_defaults(func=cmd_detect)
    submit = sub.add_parser("submit", help="normalize, queue, and detach one or more GWAS jobs")
    submit.add_argument("paths", nargs="+", help="GWAS files or directories")
    submit.add_argument("--workers", type=int, default=2)
    submit.add_argument("--run-id")
    submit.add_argument("--coloc")
    submit.set_defaults(func=cmd_submit)
    worker = sub.add_parser("worker", help="internal detached queue worker")
    worker.add_argument("--workers", type=int, default=2)
    worker.add_argument("--queue-root", default=str(QUEUE_ROOT))
    worker.set_defaults(func=lambda args: worker_loop(max(1, min(4, args.workers)), Path(args.queue_root)))
    status = sub.add_parser("status", help="show detached job status")
    status.add_argument("--run-id")
    status.set_defaults(func=cmd_status)
    report = sub.add_parser("report", help="show reproducibility report")
    report.add_argument("run_id")
    report.set_defaults(func=cmd_report)
    interface = sub.add_parser("interface", help="generate post-GWAS-compatible files")
    interface.add_argument("run_id")
    interface.add_argument("--coloc")
    interface.set_defaults(func=cmd_interface)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    for directory in (LOCAL_ROOT, JOBS_ROOT, QUEUE_ROOT):
        directory.mkdir(parents=True, exist_ok=True)
    try:
        return int(args.func(args))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
