#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
module_root="$(cd "${script_dir}/.." && pwd)"
cd "${module_root}"

output_dir="${1:-results/paper}"
config_path="${2:-config/paper_config.R}"
workers="${SLURM_CPUS_PER_TASK:-1}"
mkdir -p "${output_dir}"

status_file="${output_dir}/batch_status.txt"
start_time="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
exit_status=0
finish() {
  exit_status=$?
  end_time="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  {
    echo "start_time_utc=${start_time}"
    echo "end_time_utc=${end_time}"
    echo "workers=${workers}"
    echo "exit_status=${exit_status}"
  } > "${status_file}"
}
trap finish EXIT

{
  echo "start_time_utc=${start_time}"
  echo "r_version=$(Rscript --vanilla -e 'cat(R.version.string)')"
  Rscript --vanilla -e 'writeLines(capture.output(sessionInfo()))'
} > "${output_dir}/runtime_environment.txt"

find R config scripts tests run_all.R ../R/method.R -type f -print0 |
  sort -z |
  xargs -0 shasum > "${output_dir}/code_checksums.sha1"

for test_file in tests/test_*.R tests/smoke_test.R; do
  Rscript --vanilla "${test_file}"
done > "${output_dir}/test_suite.log" 2>&1
Rscript --vanilla scripts/01_population_check.R \
  --config="${config_path}" \
  --output-dir="${output_dir}" \
  --workers="${workers}" \
  > "${output_dir}/population.log" 2>&1
Rscript --vanilla scripts/02_run_simulation.R \
  --config="${config_path}" \
  --output-dir="${output_dir}" \
  --workers="${workers}" \
  > "${output_dir}/simulation.log" 2>&1
Rscript --vanilla scripts/03_aggregate_results.R \
  --config="${config_path}" \
  --output-dir="${output_dir}" \
  > "${output_dir}/aggregation.log" 2>&1
Rscript --vanilla scripts/04_make_markdown_tables.R \
  --config="${config_path}" \
  --output-dir="${output_dir}" \
  > "${output_dir}/tables.log" 2>&1
