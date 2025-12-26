import argparse, json
from pathlib import Path
from validators.structural import validate_structure
from validators.consistency import validate_consistency

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--strict", type=int, default=1)
    args = ap.parse_args()

    data = json.loads(Path(args.json).read_text(encoding="utf-8"))
    strict = bool(args.strict)

    s_errs = validate_structure(data)
    c_errs = validate_consistency(data, strict=strict)

    if s_errs or c_errs:
        report = {
            "structural_error_count": len(s_errs),
            "consistency_error_count": len(c_errs),
            "structural_errors": s_errs[:200],
            "consistency_errors": c_errs[:200],
        }
        print(json.dumps(report, indent=2, ensure_ascii=False))
        raise SystemExit(1)

    print("OK: validation passed")

if __name__ == "__main__":
    main()
