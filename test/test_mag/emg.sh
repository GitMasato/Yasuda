#!/bin/bash

SCRIPT_DIR=$(
  cd "$(dirname "${0}")" || exit
  pwd
)

# file and the path without extension
DEM_PATH="${SCRIPT_DIR%/*/*}/build"
DEM="emgsim"

INI_PATH="${SCRIPT_DIR}"
INI="test_mag"

INPUT_PATH="${SCRIPT_DIR}/input"
PTCL="test"
ELE=()
MAG=("B_1" "B_2" "B_3" "B_4")
OBJ=("obj_1" "obj_2")

is_STDOUT_ON=true
is_STDERR_ON=true
STDOUT_PATH="${SCRIPT_DIR}"

# process
STDOUT_ARG="${STDOUT_PATH}/stdout.txt"
STDERR_ARG="${STDOUT_PATH}/stderr.txt"
DEM_EXE="${DEM_PATH}/${DEM}"

INPUT_ARGS+=("${INI_PATH}/${INI}.inidem")
INPUT_ARGS+=("${INPUT_PATH}/${PTCL}_ptcl.csv")

for V in "${ELE[@]}"; do
  INPUT_ARGS+=("${INPUT_PATH}/${V}.eledem")
done

for V in "${MAG[@]}"; do
  INPUT_ARGS+=("${INPUT_PATH}/${V}.magdem")
done

for V in "${OBJ[@]}"; do
  INPUT_ARGS+=("${INPUT_PATH}/${V}.objdem")
done

if [[ ${is_STDOUT_ON} = true && ${is_STDERR_ON} = true ]]; then
  eval "${DEM_EXE}" "${INPUT_ARGS[@]}" >"${STDOUT_ARG}" 2>"${STDERR_ARG}"
elif [[ ${is_STDOUT_ON} = true ]]; then
  eval "${DEM_EXE}" "${INPUT_ARGS[@]}" >"${STDOUT_ARG}"
elif [[ ${is_STDERR_ON} = true ]]; then
  eval "${DEM_EXE}" "${INPUT_ARGS[@]}" 2>"${STDERR_ARG}"
else
  eval "${DEM_EXE}" "${INPUT_ARGS[@]}"
fi

if [ -e "${INI_PATH}/${INI}/OUTPUT_1/${INI}_ptcl.mgf" ]; then
  eval nkf -Lw "${INI_PATH}/${INI}/OUTPUT_1/${INI}_ptcl.mgf" \
    >"${INI_PATH}/${INI}/OUTPUT_1/${INI}_ptcl_w.mgf"
fi
