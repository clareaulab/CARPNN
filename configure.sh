#!/usr/bin/env bash
# =============================================================================
# configure.sh - one-stop path configuration for CAR-PNN
#
# Every workflow script fills in any path it has not set itself from ./config.sh
# (written by this script). A path you hardcode in a script, or export in your
# environment, always wins: config.sh only supplies defaults. Run this once after
# installing the external tools, and again whenever a path changes.
#
#   ./configure.sh                    interactive; auto-detects what it can
#   ./configure.sh --non-interactive  use env vars / --from FILE / detected values
#   ./configure.sh --from my.env      take values from a file of VAR=value lines
#   ./configure.sh --check            verify the installation (no changes)
#   ./configure.sh --set-partitions   rewrite '#SBATCH --partition=' in the scripts
#   ./configure.sh --rewrite-legacy   rewrite old hardcoded VAR="..." lines in a fork
#                                     whose scripts predate config.sh
#   ./configure.sh --no-link          do not create ~/.config/carpnn/config.sh
#   ./configure.sh --help
#
# Value precedence: detected default < previous config.sh < --from file < env var
# < what you type at the prompt.
# =============================================================================
set -uo pipefail

CARPNN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CARPNN_DIR}/config.sh"
USER_LINK="${HOME}/.config/carpnn/config.sh"
SCRIPT_DIRS=(AlphaFold Boltz LigandMPNN)   # under workflows/; the only dirs we ever rewrite

# --- variable registry --------------------------------------------------------
ORDER=(CARPNN_PYTHON BOLTZ_PATH SIF_PATH DB_FOLDER BINDCRAFT_DIR BINDCRAFT_PYTHON
       DSSP_PATH DALPHABALL_PATH LIGANDMPNN_DIR LIGANDMPNN_PYTHON CHECKPOINT_PATH
       RFDAA_AF2_PYTHON SLURM_CPU_PARTITION SLURM_GPU_PARTITION)
declare -A KIND LABEL VAL
reg() { KIND[$1]=$2; LABEL[$1]=$3; VAL[$1]=""; }
reg CARPNN_PYTHON      exe  "python in the 'carpnn' conda env (Boltz steps 03,04,06,07,09,10)"
reg BOLTZ_PATH         exe  "boltz executable (step 05)"
reg SIF_PATH           file "ColabFold singularity .sif (step 02, local MSA search)"
reg DB_FOLDER          dir  "ColabFold database directory (step 02)"
reg BINDCRAFT_DIR      dir  "BindCraft checkout - optional, only used to locate DSSP/DAlphaBall"
reg BINDCRAFT_PYTHON   exe  "python in the BindCraft env, needs pyrosetta+scipy (step 08)"
reg DSSP_PATH          exe  "DSSP executable (step 08; BindCraft ships functions/dssp)"
reg DALPHABALL_PATH    exe  "DAlphaBall.gcc (step 08; BindCraft ships functions/DAlphaBall.gcc)"
reg LIGANDMPNN_DIR     dir  "LigandMPNN checkout (SoluableMPNN scripts)"
reg LIGANDMPNN_PYTHON  exe  "python in the ligandmpnn_env conda env"
reg CHECKPOINT_PATH    file "SolubleMPNN weights, e.g. model_params/solublempnn_v_48_020.pt"
reg RFDAA_AF2_PYTHON   exe  "python in the 'mlfold' conda env (AlphaFold monomer workflow)"
reg SLURM_CPU_PARTITION str "Slurm partition(s) for CPU jobs, comma separated (optional)"
reg SLURM_GPU_PARTITION str "Slurm partition(s) for GPU jobs, comma separated (optional)"

# --- options ------------------------------------------------------------------
INTERACTIVE=true; [ -t 0 ] || INTERACTIVE=false
DO_CHECK=false; DO_PARTITIONS=false; DO_LEGACY=false; DO_LINK=true; FROM_FILE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --non-interactive) INTERACTIVE=false ;;
    --from)            FROM_FILE="${2:?--from needs a file}"; shift ;;
    --check)           DO_CHECK=true ;;
    --set-partitions)  DO_PARTITIONS=true ;;
    --rewrite-legacy)  DO_LEGACY=true ;;
    --no-link)         DO_LINK=false ;;
    -h|--help)         sed -n '3,21p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1 (see --help)" >&2; exit 2 ;;
  esac
  shift
done

say()  { printf '%s\n' "$*"; }
ok()   { printf '  [ ok ] %s\n' "$*"; }
warn() { printf '  [warn] %s\n' "$*"; }
bad()  { printf '  [FAIL] %s\n' "$*"; }

# --- detection helpers --------------------------------------------------------
CONDA_ENVS=""; CONDA_ENVS_LOADED=false
conda_prefix_of() {  # env name -> prefix dir (empty if not found)
  local n=$1 p c
  if ! $CONDA_ENVS_LOADED; then   # 'conda env list' is slow on some clusters: run it once
    CONDA_ENVS_LOADED=true
    command -v conda >/dev/null 2>&1 && CONDA_ENVS=$(conda env list 2>/dev/null)
  fi
  p=$(printf '%s\n' "$CONDA_ENVS" | awk -v n="$n" '$1==n {print $NF}' | head -n1)
  [ -n "$p" ] && [ -d "$p" ] && { echo "$p"; return; }
  for c in "${CONDA_PREFIX:+$(dirname "$CONDA_PREFIX")}" "${CONDA_PREFIX:-}/envs" "$HOME/miniconda3/envs" \
           "$HOME/anaconda3/envs" "$HOME/mambaforge/envs" "$HOME/miniforge3/envs" "$HOME/micromamba/envs"; do
    [ -n "$c" ] && [ -d "$c/$n" ] && { echo "$c/$n"; return; }
  done
}
env_bin() {          # env name, binary -> path if executable
  local p; p=$(conda_prefix_of "$1")
  [ -n "$p" ] && [ -x "$p/bin/$2" ] && echo "$p/bin/$2"
  return 0
}
first_nonempty() { local x; for x in "$@"; do [ -n "$x" ] && { echo "$x"; return; }; done; }

detect_defaults() {
  VAL[CARPNN_PYTHON]=$(first_nonempty "$(env_bin carpnn python)" \
      "$([ "${CONDA_DEFAULT_ENV:-}" = carpnn ] && echo "${CONDA_PREFIX:-}/bin/python")")
  VAL[BOLTZ_PATH]=$(first_nonempty "$(command -v boltz 2>/dev/null)" "$(env_bin boltz2 boltz)" "$(env_bin boltz boltz)")
  VAL[BINDCRAFT_PYTHON]=$(first_nonempty "$(env_bin BindCraft python)" "$(env_bin bindcraft python)")
  VAL[LIGANDMPNN_DIR]="${CARPNN_DIR}/public/LigandMPNN"
  VAL[LIGANDMPNN_PYTHON]=$(first_nonempty "$(env_bin ligandmpnn_env python)" "$(env_bin ligandmpnn python)")
  VAL[RFDAA_AF2_PYTHON]=$(env_bin mlfold python)
}

derive() {           # defaults that depend on earlier answers
  case "$1" in
    CHECKPOINT_PATH)  [ -n "${VAL[LIGANDMPNN_DIR]}" ] && echo "${VAL[LIGANDMPNN_DIR]}/model_params/solublempnn_v_48_020.pt" ;;
    DSSP_PATH)        if [ -n "${VAL[BINDCRAFT_DIR]}" ] && [ -x "${VAL[BINDCRAFT_DIR]}/functions/dssp" ]; then echo "${VAL[BINDCRAFT_DIR]}/functions/dssp"
                      else first_nonempty "$(command -v mkdssp 2>/dev/null)" "$(command -v dssp 2>/dev/null)"; fi ;;
    DALPHABALL_PATH)  [ -n "${VAL[BINDCRAFT_DIR]}" ] && [ -x "${VAL[BINDCRAFT_DIR]}/functions/DAlphaBall.gcc" ] && echo "${VAL[BINDCRAFT_DIR]}/functions/DAlphaBall.gcc" ;;
  esac
  return 0
}

load_file() {        # source FILE in a subshell and pick up the known variables
  local f=$1 v out
  [ -f "$f" ] || return 0
  for v in "${ORDER[@]}"; do
    out=$( ( unset "${ORDER[@]}" CARPNN_DIR; . "$f" >/dev/null 2>&1; printf '%s' "${!v-}" ) )
    [ -n "$out" ] && VAL[$v]=$out
  done
  return 0
}
load_env() { local v; for v in "${ORDER[@]}"; do [ -n "${!v:-}" ] && VAL[$v]=${!v}; done; return 0; }

normalize() {        # kind, value -> absolute path / resolved command
  local kind=$1 val=$2
  [ -z "$val" ] && return 0
  val=${val/#\~/$HOME}
  case "$kind" in
    str) echo "$val" ;;
    exe) if [[ "$val" != */* ]]; then command -v "$val" 2>/dev/null || echo "$val"; else realpath -m "$val"; fi ;;
    *)   realpath -m "$val" ;;
  esac
}

# --- generate config.sh -------------------------------------------------------
gather_values() {
  local v def ans hint=""
  detect_defaults
  load_file "$CONFIG_FILE"
  [ -n "$FROM_FILE" ] && { [ -f "$FROM_FILE" ] || { echo "--from file not found: $FROM_FILE" >&2; exit 2; }; load_file "$FROM_FILE"; }
  load_env
  if command -v sinfo >/dev/null 2>&1; then hint=" (available: $(sinfo -h -o '%R' 2>/dev/null | sort -u | tr '\n' ' '))"; fi
  $INTERACTIVE && say "CAR-PNN install dir: ${CARPNN_DIR}
Press Enter to accept [default]; type '-' to leave a value empty.
"
  for v in "${ORDER[@]}"; do
    [ -z "${VAL[$v]}" ] && VAL[$v]=$(derive "$v")
    if $INTERACTIVE; then
      def=${VAL[$v]}
      printf '%s%s\n' "${LABEL[$v]}" "$([[ $v == SLURM_* ]] && echo "$hint")"
      read -r -p "  $v [$def]: " ans </dev/tty || ans=""
      if [ "$ans" = "-" ]; then VAL[$v]=""; elif [ -n "$ans" ]; then VAL[$v]=$ans; fi
    fi
    VAL[$v]=$(normalize "${KIND[$v]}" "${VAL[$v]}")
  done
}

validate_var() {     # returns 0 ok, 1 warn (unset), 2 fail (set but wrong)
  local v=$1 val=${VAL[$1]}
  if [ -z "$val" ]; then
    case "$v" in
      BINDCRAFT_DIR|SLURM_*) ok "$v not set (optional)"; return 0 ;;
    esac
    warn "$v not set - ${LABEL[$v]}"; return 1
  fi
  case "${KIND[$v]}" in
    exe)  [ -x "$val" ] && { ok "$v = $val"; return 0; } ;;
    file) [ -f "$val" ] && { ok "$v = $val"; return 0; } ;;
    dir)  [ -d "$val" ] && { ok "$v = $val"; return 0; } ;;
    str)  ok "$v = $val"; return 0 ;;
  esac
  bad "$v = $val does not exist / is not executable"; return 2
}

write_config() {
  local v tmp="${CONFIG_FILE}.tmp.$$"
  {
    echo "# Generated by configure.sh on $(date -u +%Y-%m-%dT%H:%M:%SZ). Safe to edit by hand;"
    echo "# re-run ./configure.sh to regenerate. See config.example.sh for documentation."
    echo "#"
    echo "# These are DEFAULTS: a variable already set (hardcoded in a script, or exported in your"
    echo "# environment) is never overwritten by this file."
    printf '[ -n "${CARPNN_DIR:-}" ] || CARPNN_DIR=%q; export CARPNN_DIR\n' "$CARPNN_DIR"
    for v in "${ORDER[@]}"; do
      printf '\n# %s\n[ -n "${%s:-}" ] || %s=%q; export %s\n' "${LABEL[$v]}" "$v" "$v" "${VAL[$v]}" "$v"
    done
    echo
    echo "true   # keep a zero exit status when sourced"
  } > "$tmp" && mv "$tmp" "$CONFIG_FILE"
  say "Wrote $CONFIG_FILE"
}

link_user_config() {
  $DO_LINK || return 0
  local cur
  mkdir -p "$(dirname "$USER_LINK")"
  if [ -e "$USER_LINK" ] || [ -L "$USER_LINK" ]; then
    cur=$(readlink -f "$USER_LINK" 2>/dev/null || true)
    if [ "$cur" != "$(readlink -f "$CONFIG_FILE")" ]; then
      warn "$USER_LINK already points to ${cur:-a regular file}; left unchanged."
      warn "Scripts will use that other install unless you 'export CARPNN_CONFIG=$CONFIG_FILE'."
      return 0
    fi
  fi
  ln -sfn "$CONFIG_FILE" "$USER_LINK" && say "Linked $USER_LINK -> $CONFIG_FILE (lets sbatch jobs find the config)"
}

# --- in-place rewrites --------------------------------------------------------
script_files() {     # shell scripts (+ bindcraft_utils.py) we are allowed to rewrite
  local d
  for d in "${SCRIPT_DIRS[@]}"; do
    [ -d "$CARPNN_DIR/workflows/$d" ] && find "$CARPNN_DIR/workflows/$d" -type f \( -name '*.sh' -o -name 'bindcraft_utils.py' \) -not -path '*/data1/*'
  done
}
sed_escape() { printf '%s' "$1" | sed -e 's/[\\|&]/\\&/g'; }
rewrite_file() {     # file, sed-script ; returns 0 if file changed
  local f=$1 expr=$2 tmp="$1.cfgnew.$$"
  sed -e "$expr" "$f" > "$tmp"
  if cmp -s "$f" "$tmp"; then rm -f "$tmp"; return 1; fi
  cat "$tmp" > "$f"; rm -f "$tmp"; return 0
}

set_partitions() {
  local f val n=0 total=0
  say "Setting #SBATCH --partition lines (GPU scripts -> '${VAL[SLURM_GPU_PARTITION]}', others -> '${VAL[SLURM_CPU_PARTITION]}')"
  while IFS= read -r f; do
    grep -qE '^#SBATCH[[:space:]]+--partition([=[:space:]])' "$f" || continue
    total=$((total+1))
    if grep -q -- '--gres=gpu' "$f"; then val=${VAL[SLURM_GPU_PARTITION]}; else val=${VAL[SLURM_CPU_PARTITION]}; fi
    [ -z "$val" ] && { warn "no partition configured for $(basename "$f"); skipped"; continue; }
    rewrite_file "$f" "s|^#SBATCH[[:space:]]\{1,\}--partition[=[:space:]].*|#SBATCH --partition=$(sed_escape "$val")|" && n=$((n+1))
  done < <(script_files)
  say "  updated $n of $total scripts that set a partition"
}

rewrite_legacy() {
  local f v val n=0
  declare -A MAP=( [CARPNN_DIR]="$CARPNN_DIR" [CARPNN_PYTHON]="${VAL[CARPNN_PYTHON]}" [PYTHON_PATH]="${VAL[CARPNN_PYTHON]}"
    [BOLTZ_PATH]="${VAL[BOLTZ_PATH]}" [BINDCRAFT_PYTHON]="${VAL[BINDCRAFT_PYTHON]}" [DSSP_PATH]="${VAL[DSSP_PATH]}"
    [DALPHABALL_PATH]="${VAL[DALPHABALL_PATH]}" [LIGANDMPNN_DIR]="${VAL[LIGANDMPNN_DIR]}"
    [LIGANDMPNN_PYTHON]="${VAL[LIGANDMPNN_PYTHON]}" [CHECKPOINT_PATH]="${VAL[CHECKPOINT_PATH]}"
    [RFDAA_AF2_PYTHON]="${VAL[RFDAA_AF2_PYTHON]}" [SIF_PATH]="${VAL[SIF_PATH]}" [DB_FOLDER]="${VAL[DB_FOLDER]}" )
  say "Rewriting hardcoded VAR=\"...\" assignment lines in workflows/{${SCRIPT_DIRS[*]}}"
  while IFS= read -r f; do
    for v in "${!MAP[@]}"; do
      val=${MAP[$v]}
      [ -z "$val" ] && continue
      if rewrite_file "$f" "s|^\([[:space:]]*\)${v}=\"[^\"]*\"|\1${v}=\"$(sed_escape "$val")\"|"; then
        say "  ${f#$CARPNN_DIR/}: $v"; n=$((n+1))
      fi
    done
  done < <(script_files)
  say "  $n assignment(s) rewritten"
}

# --- check --------------------------------------------------------------------
run_check() {
  local v rc=0 r py
  say "Checking $CONFIG_FILE"
  [ -f "$CONFIG_FILE" ] || { bad "config.sh not found - run ./configure.sh"; return 1; }
  load_file "$CONFIG_FILE"
  for v in "${ORDER[@]}"; do validate_var "$v"; r=$?; [ $r -eq 2 ] && rc=1; done

  say "Tools"
  command -v sbatch >/dev/null 2>&1 && ok "sbatch found" || warn "sbatch not found (needed to submit any workflow job)"
  if [ -x "${VAL[CARPNN_PYTHON]}" ]; then
    py=${VAL[CARPNN_PYTHON]}
    "$py" -c "import Bio, numpy, pandas, yaml, tqdm" 2>/dev/null && ok "carpnn env imports Bio/numpy/pandas/yaml/tqdm" || { bad "carpnn env is missing Bio/numpy/pandas/yaml/tqdm (conda env create -f carpnn.yml)"; rc=1; }
  fi
  if [ -x "${VAL[BINDCRAFT_PYTHON]}" ]; then
    "${VAL[BINDCRAFT_PYTHON]}" -c "import pyrosetta, scipy, Bio, pandas" 2>/dev/null && ok "BindCraft env imports pyrosetta/scipy/Bio/pandas" || { bad "BindCraft env cannot import pyrosetta/scipy/Bio/pandas"; rc=1; }
  fi
  if [ -x "${VAL[BOLTZ_PATH]}" ]; then
    "${VAL[BOLTZ_PATH]}" predict --help >/dev/null 2>&1 && ok "boltz predict --help runs" || { bad "'boltz predict --help' failed"; rc=1; }
  fi
  [ -n "${VAL[LIGANDMPNN_DIR]}" ] && { [ -f "${VAL[LIGANDMPNN_DIR]}/run.py" ] && ok "LigandMPNN run.py found" || { bad "LigandMPNN run.py not found in ${VAL[LIGANDMPNN_DIR]}"; rc=1; }; }
  if [ -n "${VAL[SIF_PATH]}" ]; then
    command -v singularity >/dev/null 2>&1 || command -v apptainer >/dev/null 2>&1 && ok "singularity/apptainer found" || warn "neither singularity nor apptainer on PATH here (may be a cluster module)"
  fi
  [ -f "$USER_LINK" ] && ok "$USER_LINK present" || warn "$USER_LINK missing - jobs need 'export CARPNN_CONFIG=$CONFIG_FILE' (or re-run without --no-link)"

  say "Leftover hardcoded author paths"
  local pat="/data1/""lareauc" hits
  hits=$(grep -rIn --include='*.sh' --include='*.py' --include='*.md' --include='*.yml' \
         --exclude-dir=examples --exclude-dir=dev_examples --exclude-dir=public --exclude-dir=notebooks \
         --exclude-dir=.git --exclude-dir=data1 --exclude=configure.sh --exclude=config.example.sh --exclude=config.sh \
         -e "$pat" "$CARPNN_DIR" 2>/dev/null | grep -vE '^[^:]+:[0-9]+:[[:space:]]*#' | cut -c1-160)
  if [ -n "$hits" ]; then warn "found; run ./configure.sh --rewrite-legacy for scripts, edit docs by hand:"; printf '%s\n' "$hits" | sed 's/^/         /'
  else ok "none in scripts/docs (notebooks and examples are intentionally not checked)"; fi
  return $rc
}

# --- main ---------------------------------------------------------------------
if [ ! -f "$CONFIG_FILE" ] || { ! $DO_CHECK && ! $DO_PARTITIONS && ! $DO_LEGACY; }; then
  gather_values
  say ""; for v in "${ORDER[@]}"; do validate_var "$v" >/dev/null; done
  write_config
  link_user_config
  if $INTERACTIVE && [ -n "${VAL[SLURM_CPU_PARTITION]}${VAL[SLURM_GPU_PARTITION]}" ] && ! $DO_PARTITIONS; then
    read -r -p "Rewrite '#SBATCH --partition' lines in the workflow scripts now? [y/N] " a </dev/tty || a=""
    [[ $a == [yY]* ]] && DO_PARTITIONS=true
  fi
else
  load_file "$CONFIG_FILE"
fi

$DO_LEGACY     && rewrite_legacy
$DO_PARTITIONS && set_partitions

rc=0
if $DO_CHECK; then run_check; rc=$?
else say "
Next: ./configure.sh --check   (verifies every path and tool)"; fi
exit $rc
