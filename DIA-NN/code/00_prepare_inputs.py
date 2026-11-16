"""Extract DIA-NN analysis columns and validate the accompanying cleanup audit.

Input reports are read only. Configure paths with COPPRA_PACKAGE,
COPPRA_INPUT_DIR, and COPPRA_RESULTS_DIR, as in CoPPRA_config.R.
"""
import hashlib
import json
import os
from pathlib import Path

import pyarrow as pa
import pyarrow.csv as csv
import pyarrow.parquet as pq

COLUMNS = [
    "Run", "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
    "Precursor.Quantity", "Protein.Group", "Protein.Ids", "Genes",
    "Protein.Names", "Decoy", "Q.Value", "Global.Q.Value", "Global.PG.Q.Value",
    "PG.Q.Value", "Lib.Peptidoform.Q.Value",
]
INPUTS = ["report.parquet", "clean.parquet", "cleanup_audit.json", "samples.tsv",
          "contaminant_ids.txt", "clean.protein_matrix.tsv"]


def main():
    root = Path(os.environ.get("COPPRA_PACKAGE", Path(__file__).resolve().parents[1])).expanduser()
    raw = Path(os.environ.get("COPPRA_INPUT_DIR", root / "Raw_Data")).expanduser()
    results = Path(os.environ.get("COPPRA_RESULTS_DIR", root / "results")).expanduser()
    missing = [name for name in INPUTS if not (raw / name).is_file()]
    if missing:
        raise FileNotFoundError("Missing DIA-NN inputs: " + ", ".join(missing))
    audit = json.loads((raw / "cleanup_audit.json").read_text())
    report = pq.ParquetFile(raw / "report.parquet")
    clean = pq.ParquetFile(raw / "clean.parquet")
    if (raw / "report.parquet").stat().st_size != audit["source"]["bytes"]:
        raise ValueError("Source report size disagrees with cleanup_audit.json")
    if report.metadata.num_rows != audit["counts"]["input_rows"]:
        raise ValueError("Source row count disagrees with cleanup_audit.json")
    if clean.metadata.num_rows != audit["counts"]["clean_rows"]:
        raise ValueError("Clean row count disagrees with cleanup_audit.json")
    missing_columns = sorted(set(COLUMNS) - set(report.schema_arrow.names))
    if missing_columns:
        raise ValueError("Missing report columns: " + ", ".join(missing_columns))
    out = results / "Input"
    out.mkdir(parents=True, exist_ok=True)
    with pa.CompressedOutputStream(str(out / "report_analysis_columns.csv.gz"), "gzip") as sink:
        with csv.CSVWriter(sink, pa.schema([report.schema_arrow.field(c) for c in COLUMNS])) as writer:
            for batch in report.iter_batches(columns=COLUMNS, batch_size=65536):
                writer.write_batch(batch)
    manifest = []
    for name in INPUTS:
        path = raw / name
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        manifest.append({"file": name, "bytes": path.stat().st_size, "sha256": digest.hexdigest()})
    (results / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("CoPPRA: source rows", report.metadata.num_rows, "; clean rows", clean.metadata.num_rows)


if __name__ == "__main__":
    main()
