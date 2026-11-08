#!/bin/bash
#SBATCH --job-name=hardcore_vmc_dm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=norm2
#SBATCH --mem=4G

# Usage: bash /DATA/cuic/Hardcore_boson_VMC/.backups/dm_launcher_6pmf_iwv/run_double_matching.sh
# This single-file launcher snapshots the input/sources, then submits itself.
# Slurm output/error paths are set by the submission below, not fixed here.
set -euo pipefail
export TZ=Asia/Shanghai
PROJECT=/DATA/cuic/Hardcore_boson_VMC/.backups/dm_launcher_6pmf_iwv

if [[ ${1:-} != --run ]]; then
    if [[ -n ${SLURM_JOB_ID:-} ]]; then
        echo 'Start with bash run_double_matching.sh; the launcher calls sbatch automatically.' >&2
        exit 2
    fi
    if [[ $# -gt 1 ]]; then
        echo 'Usage: bash run_double_matching.sh [input.inp]'  >&2
        exit 2
    fi
    command -v sbatch >/dev/null
    # Python reads key=value or key value and strips # comments, just like C++.
    # The original input is preserved; only the copied output setting is changed.
    INPUT_FILE=${1:-$PROJECT/examples/optimize.inp}
    RUN_DIR=$(python3 - "$PROJECT" "$INPUT_FILE" <<'PY'
import datetime
import math
import pathlib
import sys

project = pathlib.Path(sys.argv[1])
source = pathlib.Path(sys.argv[2]).resolve()
text = source.read_text()
values = {}
for number, line in enumerate(text.splitlines(), 1):
    fields = line.split('#', 1)[0].replace('=', ' ').split()
    if not fields:
        continue
    if len(fields) != 2 or fields[0] in values:
        sys.exit('Invalid/duplicate input at line {}'.format(number))
    values[fields[0]] = fields[1]
required = ['Lx', 'Ly', 'Nb', 'JB', 'VB', 'K', 'chi1', 'Delta', 'delta']
for key in required:
    if key not in values:
        sys.exit('Missing parameter for run directory name: ' + key)
    if not math.isfinite(float(values[key])):
        sys.exit('Nonfinite parameter: ' + key)
lx, ly, nb = (int(values[k]) for k in ['Lx', 'Ly', 'Nb'])
if not (3 <= lx <= 128 and 3 <= ly <= 128 and lx*ly <= 4096
        and 2 <= nb and nb % 2 == 0 and nb < lx*ly):
    sys.exit('Invalid lattice/particle number for double matching')
if float(values['delta']) <= 0:
    sys.exit('delta must be strictly positive')
def label(key):
    return format(float(values[key]), '.12g')
stamp = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).strftime('%Y%m%d_%H%M%S_%f')
name = 'DM_L{}x{}_Nb{}_JB{}_VB{}_K{}_chi1{}_Delta{}_delta{}_{}'.format(
    lx, ly, nb, label('JB'), label('VB'), label('K'), label('chi1'),
    label('Delta'), label('delta'), stamp)
directory = project / 'data' / name
directory.mkdir(parents=True, exist_ok=False)
(directory / 'logs').mkdir()
(directory / 'logs' / 'optimize.original.inp').write_text(text)
lines = []
for line in text.splitlines():
    fields = line.split('#', 1)[0].replace('=', ' ').split()
    if fields and fields[0] == 'output':
        continue
    lines.append(line)
lines.extend(['', '# Managed by run_double_matching.sh; relative to this run directory.', 'output = results'])
(directory / 'logs' / 'run.inp').write_text('\n'.join(lines) + '\n')
(directory / 'logs' / 'status.txt').write_text('PREPARED\n')
print(str(directory))
PY
)
    # Shared temporary snapshot, accessible from login and compute nodes.
    # It is outside data/ and removed when the job finishes.
    mkdir -p "$PROJECT/.job_builds"
    BUILD_DIR=$(mktemp -d "$PROJECT/.job_builds/dm.XXXXXXXX")
    trap 'rc=$?; if ((rc != 0)); then printf "SUBMISSION_FAILED exit=%s\n" "$rc" > "$RUN_DIR/logs/status.txt"; case "$BUILD_DIR" in "$PROJECT"/.job_builds/dm.*) rm -rf -- "$BUILD_DIR";; esac; fi' EXIT
    cp -a "$PROJECT/src" "$PROJECT/include" "$PROJECT/tests" "$PROJECT/Makefile" "$BUILD_DIR/"
    cp "$PROJECT/run_double_matching.sh" "$RUN_DIR/logs/run_double_matching.sh"
    printf 'QUEUED\n' > "$RUN_DIR/logs/status.txt"
    SBATCH_EXTRA=()
    if [[ -n ${VMC_DEPENDENCY:-} ]]; then SBATCH_EXTRA+=(--dependency="$VMC_DEPENDENCY"); fi
    JOB_ID=$(sbatch "${SBATCH_EXTRA[@]}" --parsable --chdir="$RUN_DIR" \
        --output="$RUN_DIR/logs/slurm-%j.out" --error="$RUN_DIR/logs/slurm-%j.err" \
        "$RUN_DIR/logs/run_double_matching.sh" --run "$RUN_DIR" "$BUILD_DIR")
    printf '%s\n' "$JOB_ID" > "$RUN_DIR/logs/job_id.txt"
    printf 'Submitted job: %s\nRun directory: %s\n' "$JOB_ID" "$RUN_DIR"
    exit 0
fi

[[ $# -eq 3 && -n ${SLURM_JOB_ID:-} ]] || { echo '--run is reserved for Slurm.' >&2; exit 2; }
RUN_DIR=$2
BUILD_DIR=$3
case "$BUILD_DIR" in "$PROJECT"/.job_builds/dm.*) ;; *) echo "Invalid build directory" >&2; exit 2;; esac
cd "$RUN_DIR"
exec > >(tee -a logs/run.log) 2>&1
finish() {
    rc=$?
    rm -rf -- "$BUILD_DIR"
    if ((rc == 0)); then state=COMPLETED; else state=FAILED; fi
    printf '%s exit=%s time=%s\n' "$state" "$rc" "$(date --iso-8601=seconds)" > logs/status.txt
}
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
printf 'RUNNING job=%s time=%s\n' "$SLURM_JOB_ID" "$(date --iso-8601=seconds)" > logs/status.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
CXX=${CXX:-g++}
{
    printf 'Job ID: %s\nNode: %s\nStart: %s\n' "$SLURM_JOB_ID" "$(hostname)" "$(date --iso-8601=seconds)"
    printf 'Input snapshot: %s/logs/run.inp\n' "$RUN_DIR"
    "$CXX" --version
    printf 'Serial double matching\n'
} > logs/environment.txt
cat logs/environment.txt
echo 'Building and testing the source snapshot...'
make -C "$BUILD_DIR" CXX="$CXX" LDFLAGS=-lstdc++fs double_matching test
sha256sum logs/run.inp "$BUILD_DIR"/src/*.cpp "$BUILD_DIR"/include/*.hpp \
    "$BUILD_DIR/include/double_matching_vmc.cpp" "$BUILD_DIR/Makefile" \
    "$BUILD_DIR/hardcore_vmc_double_matching" > logs/checksums.sha256
echo 'Running double-matching VMC...'
"$BUILD_DIR/hardcore_vmc_double_matching" logs/run.inp
echo "Completed. Results: $RUN_DIR/results"
