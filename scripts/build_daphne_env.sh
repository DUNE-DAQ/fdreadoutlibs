#!/usr/bin/env bash
# Run in a clean CERN AlmaLinux 9 shell with CVMFS mounted:
# source build_daphne_env.sh /path/to/new/workarea [jobs]
# Pins the DAQ sources tested with DAPHNE-015 gateware 14f56c3 (ABI 2.0).
# The software-trigger tests used continuation off and 10-ms DAQ windows at 10 Hz.
# This builds software; it does not configure hardware or start a DAQ run.

daphne_build_env() {
    local target="${1:-$PWD/daphne-fddaq-v5.6.2}" jobs="${2:-8}"
    local repo commit extra
    local dbt_root=/cvmfs/dunedaq.opensciencegrid.org/tools/dbt/v8.14.0
    local configuration=/nfs/sw/marroyav/daphne-14f56c3/config-30462f7
    if [[ "$target" == --help || "$target" == -h ]]; then
        printf 'source build_daphne_env.sh /path/to/new/workarea [jobs]\n'
        return 0
    fi
    if [[ -n "${DBT_WORKAREA_ENV_SCRIPT_SOURCED:-}" ]]; then
        printf 'Start from a clean shell without another DAQ environment.\n' >&2
        return 1
    fi
    [[ "$jobs" =~ ^[1-9][0-9]*$ ]] || return 1
    [[ ! -e "$target" ]] || { printf 'Target already exists: %s\n' "$target" >&2; return 1; }
    [[ -r "$dbt_root/env.sh" ]] || { printf 'CVMFS DBT v8.14.0 is required.\n' >&2; return 1; }
    [[ -r "$configuration/sessions/pds-vst-session.data.xml" ]] || { printf 'Shared test configuration is unavailable: %s\n' "$configuration" >&2; return 1; }
    command -v git >/dev/null || return 1
    source "$dbt_root/env.sh" || return
    dbt-create -b stable fddaq-v5.6.2-a9-1 "$target" || return
    cd "$target" || return
    while read -r repo commit extra; do
        [[ -n "$repo" && -z "$extra" ]] || return 1
        git init -q "sourcecode/$repo" || return
        git -C "sourcecode/$repo" remote add origin "https://github.com/DUNE-DAQ/$repo.git" || return
        git -C "sourcecode/$repo" fetch --depth 1 origin "$commit" || return
        git -C "sourcecode/$repo" checkout -q --detach "$commit" || return
    done <<'SOURCES'
appmodel 3a2c8138a6e9a5828159f4de43e97148633689a5
daphnemodules c6c6ef82e9de5b9bf1379e07f3bbdb9885ea7537
dpdklibs 0711edec6b40f52a19f522f0a321128dd7836b4e
fddetdataformats 75451f42c049bf51132493075d61a1c884592425
fdreadoutlibs 9c20f875a19fe56d2597292111dd3c8f50e8d18d
fdreadoutmodules 398b621c7c448308cc82187045a43cf7a4d9546d
hermesmodules d7de1de6627dff097f282d4753508b588dad10a9
rawdatautils ab30656d57dda8daa198de3db1e44b35135a706f
SOURCES
    mkdir -p config || return
    cp -R "$configuration" config/tp-live-10ms-10Hz || return
    source env.sh || return
    dbt-build -j "$jobs" || return
    ctest --test-dir build/fdreadoutlibs --output-on-failure --no-tests=error || return
    ctest --test-dir build/dpdklibs --output-on-failure --no-tests=error || return
    source config/tp-live-10ms-10Hz/setup_db_path.sh || return
    printf 'DAPHNE environment built and loaded: %s\n' "$PWD"
}

daphne_build_env "$@"
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    exit "$?"
else
    return "$?"
fi
