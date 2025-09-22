#!/usr/bin/env bash
set -euo pipefail

# ============================
# Environment & defaults
# ============================
export LC_ALL=C
export LANG=C

EXECUTABLE="${EXECUTABLE:-/home/starlight/STARLIGHTv04/StarlightChains_v04.amd64_g77-3.4.6-r1_static.exe}"
INPUT_FILE="${INPUT_FILE:-/home/starlight/shared_directory/config_files_starlight/grid_example.in}"

# Where original spectra live (long path)
SRC="${SRC:-/home/starlight/shared_directory/config_files_starlight/spectrum}"
# Short-path where STARLIGHT will read from (symlinks live here)
DST="${DST:-/s/spec}"
# Sanitized copies go here; symlinks will point to these
CLEAN_DIR="${CLEAN_DIR:-/s/spec_clean}"
# Temp workdir
TMP_DIR="${TMP_DIR:-/tmp/starlight_run}"

VERBOSE="${VERBOSE:-1}"       # 1: print logs; 0: quiet
KEEP_TEMP="${KEEP_TEMP:-0}"   # 1: keep temp files; 0: auto-clean
PREVIEW="${PREVIEW:-0}"       # 1: show first lines of sanitized sp_0001.txt (debug)

TMP_CFG="$TMP_DIR/starlight.config"
HEADER_FILE="$TMP_DIR/header.part"
BODY_FULL="$TMP_DIR/body_full.list"   # full original body lines
BODY_FIRST="$TMP_DIR/body_first.list" # first tokens (filenames)
RENAME_MAP="$TMP_DIR/rename_map.tsv"

log() { [[ "$VERBOSE" == "1" ]] && echo -e "$*"; }

cleanup() {
  if [[ "$KEEP_TEMP" == "1" ]]; then
    log ">> Keeping temporary files (KEEP_TEMP=1)."
    return 0
  fi
  log "==> Cleaning temporary files…"
  rm -f -- "$TMP_CFG" "$HEADER_FILE" "$BODY_FULL" "$BODY_FIRST" "$RENAME_MAP" 2>/dev/null || true
  if [[ -d "$DST" ]]; then
    find "$DST" -maxdepth 1 -type l -name 'sp_*.txt' -exec rm -f {} \; 2>/dev/null || true
    rmdir "$DST" 2>/dev/null || true
  fi
  if [[ -d "$CLEAN_DIR" ]]; then
    find "$CLEAN_DIR" -maxdepth 1 -type f -name 'sp_*.txt' -delete 2>/dev/null || true
    rmdir "$CLEAN_DIR" 2>/dev/null || true
  fi
  rmdir "$TMP_DIR" 2>/dev/null || true
}
trap cleanup EXIT

# ============================
# Basic checks
# ============================
[[ -f "$EXECUTABLE" ]] || { echo "ERROR: executable not found: $EXECUTABLE" >&2; exit 1; }
[[ -f "$INPUT_FILE" ]] || { echo "ERROR: base grid not found: $INPUT_FILE" >&2; exit 1; }
[[ -d "$SRC"        ]] || { echo "ERROR: spectra directory does not exist: $SRC" >&2; exit 1; }

# ============================
# Prepare folders
# ============================
mkdir -p "$TMP_DIR" "$DST" "$CLEAN_DIR"
# Clean previous short-name symlinks & old sanitized files
find "$DST" -maxdepth 1 -type l -name 'sp_*.txt' -exec rm -f {} \; 2>/dev/null || true
find "$CLEAN_DIR" -maxdepth 1 -type f -name 'sp_*.txt' -delete 2>/dev/null || true

log "==> Preparing short-path environment for STARLIGHT"

# ============================
# Extract header ONLY up to [IsFlagSpecAvailable]
# and capture body lines separately
# ============================
: > "$HEADER_FILE"
: > "$BODY_FULL"
: > "$BODY_FIRST"

awk '
  BEGIN{in_body=0}
  {
    if (in_body==0) {
      print $0 > "'"$HEADER_FILE"'"
      if (index($0,"[IsFlagSpecAvailable]")>0) { in_body=1; next }
    } else {
      if ($0 ~ /^[[:space:]]*#/) next
      if ($0 ~ /^[[:space:]]*$/) next
      print $0 > "'"$BODY_FULL"'"
      n=split($0,a,/[[:space:]]+/)
      if (n>=1) print a[1] > "'"$BODY_FIRST"'"
    }
  }
' "$INPUT_FILE"

[[ -s "$HEADER_FILE" ]] || { echo "ERROR: could not extract header from base grid." >&2; exit 1; }
[[ -s "$BODY_FULL"  ]] || { echo "ERROR: no spectra listed in grid body (after [IsFlagSpecAvailable])." >&2; exit 1; }

# Log template preview (last non-comment, non-empty line)
if command -v tac >/dev/null 2>&1; then
  GRID_TAIL="$(tac "$INPUT_FILE" | awk 'NF && $0 !~ /^[[:space:]]*#/ {print; exit}')"
else
  GRID_TAIL="$(awk 'NF && $0 !~ /^[[:space:]]*#/ {buf=$0} END{print buf}' "$INPUT_FILE")"
fi
log ">> Using template from base grid: ${GRID_TAIL}"

# Count entries
n_body=$(wc -l < "$BODY_FIRST" | tr -d '[:space:]')
log ">> Will use exactly the spectra listed in the grid body: $n_body file(s)."

# ============================
# Sanitize each listed spectrum and make short symlinks
# ============================
: > "$RENAME_MAP"
i=1
created=0

sanitize_awk='
BEGIN{ header_done=0; }
NR==1 {
  # normalize a "lambda flux" header if present (case-insensitive)
  low=$0; for(i=1;i<=length(low);i++){ c=substr(low,i,1); if(c>="A"&&c<="Z") low=substr(low,1,i-1)"" tolower(c) "" substr(low,i+1) }
  if (low ~ /^ *lambda[[:space:]]+flux/){ print "lambda flux"; header_done=1; next }
}
{
  # trim
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0);
  if ($0=="") next;
  # split on whitespace
  n=split($0, a, /[[:space:]]+/);
  if (n<2) next;
  x=a[1]; y=a[2];
  # accept numbers like 1, 1.0, .5, -1.2e-3, +2E+4
  num="^[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?$";
  if (x ~ num && y ~ num) {
    print x, y;
  }
}
'

while IFS= read -r base; do
  base_only="${base##*/}"
  src_path="$SRC/$base_only"
  [[ -f "$src_path" ]] || { echo "ERROR: listed spectrum not found: $src_path" >&2; exit 1; }

  printf -v short "sp_%04d.txt" "$i"
  out_clean="$CLEAN_DIR/$short"

  # Sanitize: remove BOM, drop CR, replace comma with dot, keep only two numeric cols
  # If the first line is a header, normalize to exactly: "lambda flux"
  sed '1s/^\xEF\xBB\xBF//' "$src_path" \
    | tr -d '\r' \
    | sed 's/,/./g' \
    | awk "$sanitize_awk" \
    > "$out_clean"

  [[ -s "$out_clean" ]] || { echo "ERROR: sanitized output is empty for '$src_path' -> '$out_clean'." >&2; exit 1; }

  ln -sfn "$out_clean" "$DST/$short"
  printf "%s\t%s\n" "$base_only" "$short" >> "$RENAME_MAP"

  i=$((i+1)); created=$((created+1))
done < "$BODY_FIRST"

log ">> Created $created short link(s)."

# Optional preview for debugging
if [[ "$PREVIEW" == "1" && -f "$CLEAN_DIR/sp_0001.txt" ]]; then
  echo "---- Preview of sanitized sp_0001.txt ----"
  head -n 5 "$CLEAN_DIR/sp_0001.txt" || true
  echo "------------------------------------------"
fi

# ============================
# Build temporary config
# ============================
log "==> Building temporary config: $TMP_CFG"

# Start with header only and normalize CRLF
cp -f "$HEADER_FILE" "$TMP_CFG"
sed -i 's/\r$//' "$TMP_CFG"

# Point [obs_dir] to DST with trailing slash
OBS_DIR="${DST%/}/"
sed -E -i 's#^.*\[[[:space:]]*obs_dir[[:space:]]*\].*$#'"$OBS_DIR"'  [obs_dir]#' "$TMP_CFG"

# Update first line count
first_line="$(head -n1 "$TMP_CFG")"
rest_first="${first_line#* }"
echo "$created $rest_first" | cat - <(tail -n +2 "$TMP_CFG") > "$TMP_CFG.tmp" && mv "$TMP_CFG.tmp" "$TMP_CFG"

# Append rebuilt body lines: replace ONLY the first token by the short name
exec 3<"$RENAME_MAP"
while IFS= read -r orig_line && IFS=$'\t' read -r orig_base short <&3; do
  new_line=$(awk -v newfirst="$short" '
    { $1 = newfirst; out=$1; for(i=2;i<=NF;i++) out=out" "$i; print out }
  ' <<< "$orig_line")
  echo "$new_line" >> "$TMP_CFG"
done < "$BODY_FULL"
exec 3<&-

# ============================
# Run STARLIGHT
# ============================
log "==> Summary:"
log "  EXECUTABLE = $EXECUTABLE"
log "  Base grid  = $INPUT_FILE"
log "  SRC (long) = $SRC"
log "  DST (short)= $DST"
log "  TMP_CFG    = $TMP_CFG"
log

log "==> Running STARLIGHT..."
"$EXECUTABLE" < "$TMP_CFG"
ret=$?
log "==> STARLIGHT finished with exit code $ret"
exit $ret

