#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------
# Main paths (adjust if needed)
# ----------------------------------------
EXECUTABLE="/home/starlight/STARLIGHTv04/StarlightChains_v04.amd64_g77-3.4.6-r1_static.exe"
INPUT_FILE="/home/starlight/shared_directory/config_files_starlight/grid_example.in"

# ----------------------------------------
# Optional environment overrides
#   - SRC: directory containing input spectra (TXT files)
#   - DST: SHORT path inside the container for symlinks
#   - TMP_DIR: temporary directory to build the short config
#   - SPECTRA_GLOB: file pattern for spectra (default *.txt)
# ----------------------------------------
SRC="${SRC:-/home/starlight/shared_directory/config_files_starlight/spectrum}"
DST="${DST:-/s/spec}"
TMP_DIR="${TMP_DIR:-/tmp/starlight_run}"
SPECTRA_GLOB="${SPECTRA_GLOB:-*.txt}"

TMP_CFG="$TMP_DIR/starlight.config"
RENAME_MAP="$TMP_DIR/rename_map.tsv"

# ----------------------------------------
# Basic checks
# ----------------------------------------
if [[ ! -f "$EXECUTABLE" ]]; then
  echo "Error: Executable not found: $EXECUTABLE" >&2
  exit 1
fi
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: Input file (base config) not found: $INPUT_FILE" >&2
  exit 1
fi
if [[ ! -d "$SRC" ]]; then
  echo "Error: Spectra directory does not exist: $SRC" >&2
  exit 1
fi

# ----------------------------------------
# Prepare short paths and temporary config
# ----------------------------------------
echo "==> Preparing short-path environment for STARLIGHT"
mkdir -p "$DST" "$TMP_DIR"

# Clean previous symlinks in DST (do not touch real files)
find "$DST" -maxdepth 1 -type l -exec rm -f {} \; || true

echo "==> Creating short-name symlinks in: $DST"
i=1
: > "$RENAME_MAP"
shopt -s nullglob
for f in "$SRC"/$SPECTRA_GLOB; do
  base="$(basename "$f")"
  printf -v short "sp_%04d.txt" "$i"   # short names: sp_0001.txt, sp_0002.txt, ...
  ln -sfn "$f" "$DST/$short"
  printf "%s\t%s\n" "$base" "$short" >> "$RENAME_MAP"
  i=$((i+1))
done
shopt -u nullglob

if [[ ! -s "$RENAME_MAP" ]]; then
  echo "ERROR: No spectra found in '$SRC' matching pattern '$SPECTRA_GLOB'." >&2
  exit 1
fi

echo "==> Building temporary config with short paths: $TMP_CFG"
cp -f "$INPUT_FILE" "$TMP_CFG"
# Remove CRLF in case the config came from Windows
sed -i 's/\r$//' "$TMP_CFG"

# Replace long absolute prefix with short one (if config uses absolute paths under $SRC/)
sed -i "s#${SRC%/}/#${DST%/}/#g" "$TMP_CFG"

# Replace long basenames with short ones according to the map (if config lists only filenames)
while IFS=$'\t' read -r long short; do
  sed -i "s#$long#$short#g" "$TMP_CFG"
done < "$RENAME_MAP"

echo "==> Summary:"
echo "  EXECUTABLE = $EXECUTABLE"
echo "  Base config = $INPUT_FILE"
echo "  SRC (long)  = $SRC"
echo "  DST (short) = $DST"
echo "  TMP_CFG     = $TMP_CFG"
echo

# ----------------------------------------
# Run STARLIGHT with the temporary config
# ----------------------------------------
echo "==> Running STARLIGHT..."
"$EXECUTABLE" < "$TMP_CFG"
ret=$?
echo "==> STARLIGHT finished with exit code $ret"
exit $ret
