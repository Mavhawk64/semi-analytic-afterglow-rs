import os
import re
from pathlib import Path

CWD = os.path.dirname(os.path.abspath(__file__))

# -------- CONFIG --------
SCRIPT_DIR = Path(__file__).resolve().parent
FORTRAN_ROOT = SCRIPT_DIR.parent / "semi-analytic-afterglow"
RUST_ROOT = SCRIPT_DIR.parent / "semi-analytic-afterglow-rs" / "src"

# Optional overrides (only if needed)
OVERRIDES = {}
EXCLUDE_DIRS = {".venv", ".git", "target", "__pycache__"}
EXCLUDE_FILES = ["bw_euler_stepper-conflicting.f90"]


# -------- HELPERS --------
def camel_to_snake(name: str) -> str:
    name = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", name)
    name = re.sub(r"_+", "_", name)

    return name.lower()


def normalize(name: str) -> str:
    if name in OVERRIDES:
        return OVERRIDES[name]
    return camel_to_snake(name)


def get_fortran_files(root: Path):
    files = []
    for path in root.rglob("*.f90"):
        if (
            any(part in EXCLUDE_DIRS for part in path.parts)
            or path.name in EXCLUDE_FILES
        ):
            continue
        files.append(path)
    return files


def get_rust_files(root: Path):
    files = []
    for path in root.rglob("*.rs"):
        if path.name == "mod.rs" or path.name == "lib.rs" or path.name == "main.rs":
            continue
        files.append(path)
    return files


# -------- MAIN --------
def main():
    f_files = get_fortran_files(FORTRAN_ROOT)
    r_files = get_rust_files(RUST_ROOT)
    missing_files = []

    f_map = {}
    for f in f_files:
        base = f.stem
        norm = normalize(base)
        f_map[norm] = f

    r_map = {}
    for r in r_files:
        r_map[r.stem] = r

    matched = []
    missing = []
    extra = []

    for norm, f_path in f_map.items():
        if norm in r_map:
            matched.append((f_path, r_map[norm]))
        else:
            missing.append((f_path, norm))

    for r_name, r_path in r_map.items():
        if r_name not in f_map:
            extra.append(r_path)

    # -------- OUTPUT --------
    print("\n=== MATCHED ===")
    for f, r in matched:
        print(f"✔ {f} -> {r}")

    print("\n=== MISSING (need Rust file) ===")
    for f, expected in missing:
        expected_file = f"{RUST_ROOT}/{str(f.parent.relative_to(FORTRAN_ROOT)) + '/' if f.parent.relative_to(FORTRAN_ROOT) != Path('.') else ''}{expected}.rs"
        missing_files.append(expected_file)
        print(f"✘ {f} -> expected {expected_file}")

    print("\n=== EXTRA (no Fortran source) ===")
    for r in extra:
        print(f"⚠ {r}")

    print("\nSummary:")
    print(f"Matched: {len(matched)}")
    print(f"Missing: {len(missing)}")
    print(f"Extra:   {len(extra)}")
    if missing_files:
        print("\n\033[93mCreate the following missing Rust files?\033[0m")
        for mf in missing_files:
            print(f"  - {mf}")
        x = input("Create missing files? (y/n) ")
        if x.lower() == "y":
            for mf in missing_files:
                mf_path = Path(mf)
                mf_path.parent.mkdir(parents=True, exist_ok=True)
                mf_path.touch()
            print("Missing files created.")
        else:
            print("No files created.")


if __name__ == "__main__":
    main()
