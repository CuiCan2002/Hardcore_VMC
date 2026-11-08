#!/bin/bash
#SBATCH --job-name=hardcore_vmc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=norm2
#SBATCH --time=24:00:00
#SBATCH --mem=4G

# Usage: bash /DATA/cuic/Hardcore_boson_VMC/run_optimize.sh
# This single-file launcher snapshots the input/sources, then submits itself.
# Slurm output/error paths are set by the submission below, not fixed here.
set -euo pipefail
export TZ=Asia/Shanghai
PROJECT=/DATA/cuic/Hardcore_boson_VMC

if [[ ${1:-} != --run ]]; then
    if [[ -n ${SLURM_JOB_ID:-} ]]; then
        echo 'Start with bash run_optimize.sh; the launcher calls sbatch automatically.' >&2
        exit 2
    fi
    if [[ $# -ne 0 ]]; then
        echo 'Usage: bash run_optimize.sh (reads examples/optimize.inp)' >&2
        exit 2
    fi
    command -v sbatch >/dev/null
    # Python reads key=value or key value and strips # comments, just like C++.
    # The original input is preserved; only the copied output setting is changed.
    RUN_DIR=$(python3 - "$PROJECT" <<'PY'
import datetime
import math
import pathlib
import sys

project = pathlib.Path(sys.argv[1])
source = project / 'examples' / 'optimize.inp'
text = source.read_text()
values = {}
for number, line in enumerate(text.splitlines(), 1):
    fields = line.split('#', 1)[0].replace('=', ' ').split()
    if not fields:
        continue
    if len(fields) != 2 or fields[0] in values:
        sys.exit('Invalid/duplicate input at line {}'.format(number))
    values[fields[0]] = fields[1]
required = ['Lx', 'Ly', 'Nb', 'JB', 'VB', 'K', 'chi1', 'Delta', 'lambda_eff']
for key in required:
    if key not in values:
        sys.exit('Missing parameter for run directory name: ' + key)
    if not math.isfinite(float(values[key])):
        sys.exit('Nonfinite parameter: ' + key)
lx, ly, nb = (int(values[k]) for k in ['Lx', 'Ly', 'Nb'])
if not (3 <= lx <= 128 and 3 <= ly <= 128 and lx*ly <= 4096
        and 2 <= nb <= 20 and nb % 2 == 0 and nb < lx*ly):
    sys.exit('Invalid lattice/particle number for this exact-Hafnian implementation')
def label(key):
    return format(float(values[key]), '.12g')
stamp = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).strftime('%Y%m%d_%H%M%S_%f')
name = 'L{}x{}_Nb{}_JB{}_VB{}_K{}_chi1{}_Delta{}_le{}_{}'.format(
    lx, ly, nb, label('JB'), label('VB'), label('K'), label('chi1'),
    label('Delta'), label('lambda_eff'), stamp)
directory = project / 'data' / name
directory.mkdir(parents=True, exist_ok=False)
(directory / 'optimize.original.inp').write_text(text)
lines = []
for line in text.splitlines():
    fields = line.split('#', 1)[0].replace('=', ' ').split()
    if fields and fields[0] == 'output':
        continue
    lines.append(line)
lines.extend(['', '# Managed by run_optimize.sh; relative to this run directory.', 'output = results'])
(directory / 'run.inp').write_text('\n'.join(lines) + '\n')
(directory / 'status.txt').write_text('PREPARED\n')
print(str(directory))
PY
)
    # Compile a per-job source snapshot, so simultaneous jobs cannot overwrite binaries.
    trap 'rc=$?; if ((rc != 0)); then printf "SUBMISSION_FAILED exit=%s\n" "$rc" > "$RUN_DIR/status.txt"; fi' EXIT
    mkdir "$RUN_DIR/source"
    cp -a "$PROJECT/src" "$PROJECT/include" "$PROJECT/tests" "$PROJECT/Makefile" "$RUN_DIR/source/"
    cp "$PROJECT/run_optimize.sh" "$RUN_DIR/run_optimize.sh"
    printf 'QUEUED\n' > "$RUN_DIR/status.txt"
    JOB_ID=$(sbatch --parsable --chdir="$RUN_DIR" \
        --output="$RUN_DIR/slurm-%j.out" --error="$RUN_DIR/slurm-%j.err" \
        "$RUN_DIR/run_optimize.sh" --run "$RUN_DIR")
    printf '%s\n' "$JOB_ID" > "$RUN_DIR/job_id.txt"
    printf 'Submitted job: %s\nRun directory: %s\n' "$JOB_ID" "$RUN_DIR"
    exit 0
fi

[[ $# -eq 2 && -n ${SLURM_JOB_ID:-} ]] || { echo '--run is reserved for Slurm.' >&2; exit 2; }
RUN_DIR=$2
cd "$RUN_DIR"
exec > >(tee -a run.log) 2>&1
finish() {
    rc=$?
    if ((rc == 0)); then state=COMPLETED; else state=FAILED; fi
    printf '%s exit=%s time=%s\n' "$state" "$rc" "$(date --iso-8601=seconds)" > status.txt
}
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
printf 'RUNNING job=%s time=%s\n' "$SLURM_JOB_ID" "$(date --iso-8601=seconds)" > status.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
CXX=${CXX:-g++}
{
    printf 'Job ID: %s\nNode: %s\nStart: %s\n' "$SLURM_JOB_ID" "$(hostname)" "$(date --iso-8601=seconds)"
    printf 'Input snapshot: %s/run.inp\n' "$RUN_DIR"
    "$CXX" --version
} > environment.txt
cat environment.txt
echo 'Building and testing the source snapshot...'

# GCC 7/8 require -lstdc++fs. Probe linking instead of relying on version strings.
cat > source/filesystem_probe.cpp <<'CPP'
#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif
int main() { return fs::exists(".") ? 0 : 1; }
CPP
FS_FLAGS=
if ! "$CXX" -std=c++17 source/filesystem_probe.cpp -o source/filesystem_probe > filesystem_probe.log 2>&1; then
    FS_FLAGS=-lstdc++fs
    "$CXX" -std=c++17 source/filesystem_probe.cpp -o source/filesystem_probe "$FS_FLAGS" >> filesystem_probe.log 2>&1
fi
printf 'Filesystem linker flags: %s\n' "$FS_FLAGS" >> environment.txt
make -C source CXX="$CXX" LDFLAGS="$FS_FLAGS" all test
sha256sum run.inp source/src/*.cpp source/include/*.hpp source/hardcore_vmc > checksums.sha256
echo 'Running VMC...'
srun --ntasks=1 --cpus-per-task=1 ./source/hardcore_vmc run.inp
echo "Completed. Results: $RUN_DIR/results"
