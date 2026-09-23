#!/usr/bin/env bash
set -Eeuo pipefail

# Task 001 resumable public-matrix downloader.
# Run this on the compute server, not on the control laptop:
#   bash scripts/shell/task001_download_public_matrices.sh
# Set DOWNLOAD_NORMALIZED_ONLY=1 to fetch the GSE290298 normalized-only file.

SERVER_ROOT="${SERVER_ROOT:-/data/lf_data/HCC_Peritumoral_Neutrophil_scRNA_Atlas}"
RAW_ROOT="${SERVER_ROOT}/raw_data"
LOG_ROOT="${SERVER_ROOT}/logs"
CHECKSUM_ROOT="${SERVER_ROOT}/checksums"
RUN_LOG="${LOG_ROOT}/task001_download_public_matrices.log"
DOWNLOAD_NORMALIZED_ONLY="${DOWNLOAD_NORMALIZED_ONLY:-0}"
DOWNLOAD_PARALLELISM="${DOWNLOAD_PARALLELISM:-12}"
DOWNLOAD_CHUNK_BYTES="${DOWNLOAD_CHUNK_BYTES:-33554432}"

mkdir -p "${RAW_ROOT}" "${LOG_ROOT}" "${CHECKSUM_ROOT}"
exec > >(tee -a "${RUN_LOG}") 2>&1

timestamp() { date -u '+%Y-%m-%dT%H:%M:%SZ'; }

file_size() {
  if stat -c '%s' "$1" >/dev/null 2>&1; then
    stat -c '%s' "$1"
  else
    stat -f '%z' "$1"
  fi
}

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1"
  else
    shasum -a 256 "$1"
  fi
}

download_one() {
  local dataset="$1"
  local filename="$2"
  local expected_bytes="$3"
  local url="$4"
  local outdir="${RAW_ROOT}/${dataset}"
  local outfile="${outdir}/${filename}"
  local logfile="${LOG_ROOT}/task001_${dataset}.log"
  local parts_dir="${outfile}.parts"

  mkdir -p "${outdir}"
  {
    echo "[$(timestamp)] START dataset=${dataset} file=${filename} expected_bytes=${expected_bytes} url=${url}"
    local existing_bytes=0
    if [[ -f "${outfile}" ]]; then
      existing_bytes="$(file_size "${outfile}")"
    fi
    if (( existing_bytes > expected_bytes )); then
      echo "[$(timestamp)] DOWNLOAD_FAILED path=${outfile} reason=local_file_larger_than_expected" >&2
      return 1
    fi
    if (( existing_bytes == expected_bytes )); then
      echo "[$(timestamp)] RESUME_OK path=${outfile} bytes=${existing_bytes}"
    else
      mkdir -p "${parts_dir}"
      local start="${existing_bytes}"
      local part_index=0
      local running=0
      local part_path part_end part_bytes
      while (( start < expected_bytes )); do
        part_end=$(( start + DOWNLOAD_CHUNK_BYTES - 1 ))
        if (( part_end >= expected_bytes )); then
          part_end=$(( expected_bytes - 1 ))
        fi
        part_path="${parts_dir}/part_$(printf '%08d' "${part_index}")"
        part_bytes=$(( part_end - start + 1 ))
        if [[ -f "${part_path}" ]] && [[ "$(file_size "${part_path}")" -eq "${part_bytes}" ]]; then
          echo "[$(timestamp)] RANGE_REUSE file=${filename} start=${start} end=${part_end} bytes=${part_bytes}"
        else
          rm -f "${part_path}"
          curl --fail --silent --show-error --location --retry 5 --retry-delay 10 --retry-all-errors \
            --connect-timeout 30 --range "${start}-${part_end}" --output "${part_path}" "${url}" &
          running=$(( running + 1 ))
        fi
        if (( running >= DOWNLOAD_PARALLELISM )); then
          wait
          running=0
        fi
        start=$(( part_end + 1 ))
        part_index=$(( part_index + 1 ))
      done
      if (( running > 0 )); then
        wait
      fi
      : > "${outfile}.append"
      for part_path in "${parts_dir}"/part_*; do
        [[ -f "${part_path}" ]] || continue
        cat "${part_path}" >> "${outfile}.append"
      done
      cat "${outfile}.append" >> "${outfile}"
      rm -f "${outfile}.append"
    fi
    if [[ "$(file_size "${outfile}")" -eq "${expected_bytes}" ]]; then
      echo "[$(timestamp)] DOWNLOAD_OK path=${outfile} bytes=$(wc -c < "${outfile}")"
      sha256_file "${outfile}" | tee -a "${CHECKSUM_ROOT}/task001_download_checksums.sha256"
      case "${filename}" in
        *.gz) gzip -t "${outfile}" && echo "[$(timestamp)] GZIP_OK path=${outfile}" ;;
        *.tar) tar -tf "${outfile}" >/dev/null && echo "[$(timestamp)] TAR_OK path=${outfile}" ;;
        *.csv.gz|*.txt.gz|*.mtx.gz|*.tsv.gz) gzip -t "${outfile}" && echo "[$(timestamp)] MATRIX_STREAM_OK path=${outfile}" ;;
      esac
    else
      echo "[$(timestamp)] DOWNLOAD_FAILED dataset=${dataset} file=${filename} reason=assembled_size_mismatch" >&2
      return 1
    fi
    echo "[$(timestamp)] END dataset=${dataset} file=${filename}"
  } 2>&1 | tee -a "${logfile}"
}

download_one GSE282701 GSE282701_RAW.tar 856524800 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282701/suppl/GSE282701_RAW.tar

download_one GSE242889 GSE242889_RAW.tar 1376245760 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242889/suppl/GSE242889_RAW.tar

download_one GSE326201 GSE326201_RAW.tar 244039680 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE326nnn/GSE326201/suppl/GSE326201_RAW.tar

download_one GSE149614 GSE149614_HCC.scRNAseq.S71915.count.txt.gz 165349783 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz
download_one GSE149614 GSE149614_HCC.metadata.updated.txt.gz 489882 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz

download_one GSE299340 GSE299340_RAW.tar 843857920 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE299nnn/GSE299340/suppl/GSE299340_RAW.tar

download_one GSE202642 GSE202642_matrix.mtx.gz 698853827 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202642/suppl/GSE202642_matrix.mtx.gz
download_one GSE202642 GSE202642_features.tsv.gz 333437 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202642/suppl/GSE202642_features.tsv.gz
download_one GSE202642 GSE202642_barcodes.tsv.gz 583228 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202642/suppl/GSE202642_barcodes.tsv.gz

if [[ "${DOWNLOAD_NORMALIZED_ONLY}" == "1" ]]; then
  download_one GSE290298 GSE290298_All_processed_data.csv.gz 923544964 \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE290nnn/GSE290298/suppl/GSE290298_All_processed_data.csv.gz
else
  echo "[$(timestamp)] SKIP GSE290298 normalized-only matrix; set DOWNLOAD_NORMALIZED_ONLY=1 for optional download"
fi

echo "[$(timestamp)] Task 001 matrix download phase completed"
