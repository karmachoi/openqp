#!/bin/bash
#SBATCH --nodes=1

set -euo pipefail

campaign_root=${OQP_URACIL_CAMPAIGN_ROOT:?set OQP_URACIL_CAMPAIGN_ROOT}
phase=${OQP_URACIL_PHASE:?set OQP_URACIL_PHASE}
threads=${OQP_OMP_THREADS:?set OQP_OMP_THREADS}
pack_base=${OQP_PACK_BASE:?set OQP_PACK_BASE}
pack_size=${OQP_PACK_SIZE:?set OQP_PACK_SIZE}
last_index=${OQP_LAST_INDEX:-49}
pack_index=${SLURM_ARRAY_TASK_ID:?run as a Slurm array}
worker="$campaign_root/scripts/run_slurm_namd_array.sbatch"

indices=()
if [ -n "${OQP_PACK_INDEX_LIST:-}" ]; then
  normalized_index_list=${OQP_PACK_INDEX_LIST//:/,}
  IFS=, read -ra explicit_indices <<< "$normalized_index_list"
  first_position=$((pack_index * pack_size))
  for ((slot = 0; slot < pack_size; slot++)); do
    position=$((first_position + slot))
    if [ "$position" -lt "${#explicit_indices[@]}" ]; then
      indices+=("${explicit_indices[$position]}")
    fi
  done
else
  start=$((pack_base + pack_index * pack_size))
  stop=$((start + pack_size - 1))
  if [ "$stop" -gt "$last_index" ]; then
    stop=$last_index
  fi
  if [ "$start" -le "$last_index" ]; then
    mapfile -t indices < <(seq "$start" "$stop")
  fi
fi

test -x "$worker" || test -f "$worker"
if [ "${#indices[@]}" -eq 0 ]; then
  exit 0
fi

pids=()
launch_mode=${OQP_PACK_LAUNCH:-srun}
if [ "$launch_mode" = direct ]; then
  declare -A allowed_cpu=()
  allowed_list=$(awk '/Cpus_allowed_list/ {print $2}' /proc/self/status)
  IFS=, read -ra allowed_parts <<< "$allowed_list"
  for part in "${allowed_parts[@]}"; do
    if [[ "$part" == *-* ]]; then
      for cpu in $(seq "${part%-*}" "${part#*-}"); do
        allowed_cpu[$cpu]=1
      done
    else
      allowed_cpu[$part]=1
    fi
  done
  declare -A seen_core=()
  physical_cpus=()
  while IFS=, read -r cpu core; do
    [[ "$cpu" == \#* ]] && continue
    if [[ -n ${allowed_cpu[$cpu]:-} && -z ${seen_core[$core]:-} ]]; then
      physical_cpus+=("$cpu")
      seen_core[$core]=1
    fi
  done < <(lscpu -p=CPU,CORE)
  printf 'packed_preflight host=%s job=%s array=%s allowed=' "$(hostname)" "$SLURM_JOB_ID" "$pack_index"
  taskset -pc $$
  printf 'physical_cpu_count=%s required=%s cpus=%s indices=%s\n' \
    "${#physical_cpus[@]}" "$((${#indices[@]} * threads))" "${physical_cpus[*]}" "${indices[*]}"
  test "${#physical_cpus[@]}" -ge $((${#indices[@]} * threads))
  if [ "${OQP_PACK_PREFLIGHT:-0}" = 1 ]; then
    exit 0
  fi
fi
slot=0
for index in "${indices[@]}"; do
  if [ "$launch_mode" = direct ]; then
    first=$((slot * threads))
    cpu_list=$(IFS=,; echo "${physical_cpus[*]:first:threads}")
    if [ "${OQP_WORKER_TRACE:-0}" = 1 ]; then
      env SLURM_ARRAY_TASK_ID="$index" SLURM_CPUS_PER_TASK="$threads" \
        OQP_URACIL_PHASE="$phase" OQP_OMP_THREADS="$threads" \
        taskset -c "$cpu_list" bash -x "$worker" &
    else
      env SLURM_ARRAY_TASK_ID="$index" SLURM_CPUS_PER_TASK="$threads" \
        OQP_URACIL_PHASE="$phase" OQP_OMP_THREADS="$threads" \
        taskset -c "$cpu_list" bash "$worker" &
    fi
  else
    srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task="$threads" \
      --cpu-bind=cores \
      env SLURM_ARRAY_TASK_ID="$index" SLURM_CPUS_PER_TASK="$threads" \
        OQP_URACIL_PHASE="$phase" OQP_OMP_THREADS="$threads" \
        bash "$worker" &
  fi
  pids+=("$!")
  slot=$((slot + 1))
done

status=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    status=1
  fi
done
exit "$status"
