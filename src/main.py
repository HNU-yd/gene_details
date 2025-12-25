
import argparse
from pathlib import Path
from species_runner import run_species

def main(data_root, out_dir, log_dir):
    for species_dir in Path(data_root).iterdir():
        if species_dir.is_dir():
            run_species(species_dir, out_dir, log_dir)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--log-dir", required=True)
    args = ap.parse_args()
    main(args.data_root, args.out, args.log_dir)
