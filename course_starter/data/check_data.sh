#!/usr/bin/env bash
# check_data.sh -- verify the course data files against data/checksums.md5.
#
# Usage (from anywhere):   bash data/check_data.sh
# Exit code 0 = all files present and intact, 1 = something is missing or corrupt.
#
# Works on Linux (md5sum) and macOS (md5). On Windows use: Rscript data/check_data.R

set -u
cd "$(dirname "$0")" || exit 1

if [ ! -f checksums.md5 ]; then
  echo "checksums.md5 not found next to this script." >&2
  exit 1
fi

# Pick a tool that prints just the md5 hash of one file.
if command -v md5sum >/dev/null 2>&1; then
  hash_of() { md5sum "$1" | cut -d ' ' -f 1; }
elif command -v md5 >/dev/null 2>&1; then
  hash_of() { md5 -q "$1"; }
else
  echo "Neither md5sum nor md5 is available. Use: Rscript data/check_data.R" >&2
  exit 1
fi

ok=0; bad=0; missing=0
while read -r expected file; do
  [ -z "${file:-}" ] && continue
  if [ ! -f "$file" ]; then
    echo "MISSING  $file"
    missing=$((missing + 1))
  elif [ "$(hash_of "$file")" = "$expected" ]; then
    echo "OK       $file"
    ok=$((ok + 1))
  else
    echo "CORRUPT  $file   (checksum differs: re-download or re-unzip)"
    bad=$((bad + 1))
  fi
done < checksums.md5

echo
echo "$ok OK, $missing missing, $bad corrupt"
if [ $((missing + bad)) -gt 0 ]; then
  echo "Data check FAILED. Extract the course data zip into data/ (so that data/raw/, data/derived/, data/reference/ exist) and run again."
  exit 1
fi
echo "Data check passed. Next: open analyses/01_explore_data.qmd"
