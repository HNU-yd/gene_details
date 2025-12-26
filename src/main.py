import argparse, json
from pathlib import Path
from build_tree import build_species_tree
from validators.structural import validate_structure
from validators.consistency import validate_consistency

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--strict", type=int, default=1)
    args = ap.parse_args()

    data_root = Path(args.data_root)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    strict = bool(args.strict)

    for species_dir in sorted([p for p in data_root.iterdir() if p.is_dir()]):
        tree = build_species_tree(species_dir, strict=strict).to_dict()

        s_errs = validate_structure(tree)
        c_errs = validate_consistency(tree, strict=strict)
        if s_errs or c_errs:
            msg = {
                "species": species_dir.name,
                "structural_error_count": len(s_errs),
                "consistency_error_count": len(c_errs),
                "structural_errors": s_errs[:200],
                "consistency_errors": c_errs[:200],
            }
            raise RuntimeError("Validation failed:\\n" + json.dumps(msg, indent=2, ensure_ascii=False))

        (out_dir / f"{species_dir.name}.json").write_text(json.dumps(tree, indent=2, ensure_ascii=False), encoding="utf-8")

if __name__ == "__main__":
    main()
