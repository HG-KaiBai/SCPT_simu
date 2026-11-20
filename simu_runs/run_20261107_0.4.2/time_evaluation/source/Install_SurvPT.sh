#!/usr/bin/env bash
set -euo pipefail
simu_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$simu_dir/Environment.sh"
library_dir=/home/y679zhao/R/empt-survpt-0.4.2
archive="$simu_dir/../SurvPT_0.4.2.tar.gz"
if [[ ! -f $archive ]]; then archive="$simu_dir/../packages/SurvPT_0.4.2.tar.gz"; fi
mkdir -p "$library_dir"
archive_md5=$(md5sum "$archive" | awk '{print $1}')
R CMD INSTALL --library="$library_dir" "$archive"
[[ $(md5sum "$archive" | awk '{print $1}') == "$archive_md5" ]] || { echo 'Archive changed during installation' >&2; exit 1; }
printf '%s\n' "$archive_md5" > "$library_dir/SurvPT/EMPT_ARCHIVE_MD5"
Rscript -e 'stopifnot(as.character(packageVersion("SurvPT")) == "0.4.2"); stopifnot(is.function(SurvPT::data_mixture)); cat(find.package("SurvPT"), "\n")'
