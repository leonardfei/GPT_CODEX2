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

mkdir -p "${RAW_ROOT}" "${LOG_ROOT}" "${CHECKSUM_ROOT}"
exec > >(tee -a "${RUN_LOG}") 2>&1

timestamp() { date -u '+%Y-%m-%dT%H:%M:%SZ'; }

download_one() {
  local dataset="$1"
  local filename="$2"
  local url="$3"
  local outdir="${RAW_ROOT}/${dataset}"
  local outfile="${outdir}/${filename}"
  local logfile="${LOG_ROOT}/task001_${dataset}.log"

  mkdir -p "${outdir}"
  {
    echo "[$(timestamp)] START dataset=${dataset} file=${filename} url=${url}"
    if curl --fail --location --retry 5 --retry-delay 10 --retry-all-errors \
      --connect-timeout 30 --continue-at - --output "${outfile}" "${url}"; then
      echo "[$(timestamp)] DOWNLOAD_OK path=${outfile} bytes=$(wc -c < "${outfile}")"
      sha256sum "${outfile}" | tee -a "${CHECKSUM_ROOT}/task001_download_checksums.sha256"
      case "${filename}" in
        *.gz) gzip -t "${outfile}" && echo "[$(timestamp)] GZIP_OK path=${outfile}" ;;
        *.tar) tar -tf "${outfile}" >/dev/null && echo "[$(timestamp)] TAR_OK path=${outfile}" ;;
        *.csv.gz|*.txt.gz|*.mtx.gz|*.tsv.gz) gzip -t "${outfile}" && echo "[$(timestamp)] MATRIX_STREAM_OK path=${outfile}" ;;
      esac
    else
      echo "[$(timestamp)] DOWNLOAD_FAILED dataset=${dataset} file=${filename}" >&2
      return 1
    fi
    echo "[$(timestamp)] END dataset=${dataset} file=${filename}"
  } 2>&1 | tee -a "${logfile}"
}

download_one GSE282701 GSE282701_RAW.tar \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282701/suppl/GSE282701_RAW.tar

download_one GSE242889 GSE242889_RAW.tar \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242889/suppl/GSE242889_RAW.tar

download_one GSE326201 GSE326201_RAW.tar \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE326nnn/GSE326201/suppl/GSE326201_RAW.tar

download_one GSE149614 GSE149614_HCC.scRNAseq.S71915.count.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz
download_one GSE149614 GSE149614_HCC.metadata.updated.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz

download_one GSE299340 GSE299340_RAW.tar \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE299nnn/GSE299340/suppl/GSE299340_RAW.tar

download_one GSE202642 GSE202642_matrix.mtx.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202642/suppl/GSE202642_matrix.mtx.gz
download_one GSE202642 GSE202642_features.tsv.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202642/suppl/GSE202642_features.tsv.gz
download_one GSE202642 GSE202642_barcodes.tsv.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202642/suppl/GSE202642_barcodes.tsv.gz

if [[ "${DOWNLOAD_NORMALIZED_ONLY}" == "1" ]]; then
  download_one GSE290298 GSE290298_All_processed_data.csv.gz \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE290nnn/GSE290298/suppl/GSE290298_All_processed_data.csv.gz
else
  echo "[$(timestamp)] SKIP GSE290298 normalized-only matrix; set DOWNLOAD_NORMALIZED_ONLY=1 for optional download"
fi

echo "[$(timestamp)] Task 001 matrix download phase completed"
