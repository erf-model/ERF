#!/usr/bin/env bash
# Linux only: Docker and GNU coreutils are available on the Ubuntu CI runners.
# Invalid inputs fail; cache transport failures warn and leave CI free to continue.
set -euo pipefail

error() { echo "::error::$*" >&2; exit 2; }
warning() {
  if [[ -f ${log:-} ]]; then tail -n 20 "$log" >&2; fi
  echo "::warning::$*" >&2
}

case "$CACHE_MODE" in restore|save) ;; *) error 'mode must be restore or save' ;; esac
[[ "$RUNNER_OS" == Linux ]] || error 'ccache-trunk requires a Linux runner'
[[ "$CACHE_MODE" != save || -n "$GHCR_TOKEN" ]] || error 'token is required when saving'
owner=$(printf '%s' "$GITHUB_REPOSITORY_OWNER" | tr '[:upper:]' '[:lower:]')
tag=$(printf '%s' "$CACHE_KEY" | sed -E 's/[^A-Za-z0-9_.-]+/-/g' | tr '[:upper:]' '[:lower:]')
[[ "$tag" =~ ^[a-z0-9_][a-z0-9_.-]{0,127}$ ]] || error 'key must produce a valid GHCR tag'
ref="ghcr.io/$owner/erf-ccache:$tag"
case "$CACHE_PATH" in
  \~) path=$HOME ;;
  \~/*) path="$HOME/${CACHE_PATH:2}" ;;
  '') error 'path must name a dedicated ccache directory' ;;
  *) path=$CACHE_PATH ;;
esac
path=$(realpath -m -- "$path")
[[ "$path" != / && "$path" != "$HOME" ]] || error 'path must name a dedicated ccache directory'
[[ ! -e "$path" || -d "$path" ]] || error 'cache path is not a directory'
if [[ "$CACHE_MODE" == restore ]]; then
  echo 'hit=false' >> "$GITHUB_OUTPUT"
elif [[ ! -d "$path" || -z $(find "$path" -type f -print -quit) ]]; then
  echo "no ccache files in $path: skip save"
  exit 0
fi

if ! scratch=$(mktemp -d "$RUNNER_TEMP/erf-ccache.XXXXXX"); then
  warning 'could not prepare GHCR cache transfer; continuing without it'
  exit 0
fi
log="$scratch/transfer.log"
container=''
stage=''
backup=''
cleanup() {
  if [[ -n "$container" ]]; then
    timeout --kill-after=5s 10s docker rm "$container" >/dev/null 2>&1 || true
  fi
  # Restore the original directory if interrupted between the two renames below.
  if [[ -n "$backup" && -e "$backup" ]]; then
    if [[ ! -e "$path" ]]; then
      mv -- "$backup" "$path" || echo "::warning::original cache remains at $backup" >&2
    else
      rm -rf -- "$backup" || true
    fi
  fi
  if [[ -n "$stage" ]]; then rm -rf -- "$stage" || true; fi
  rm -rf -- "$scratch" || true
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
# Keep credentials out of the runner's shared Docker config and remove them on exit.
export DOCKER_CONFIG="$scratch/docker"
if ! mkdir -p "$DOCKER_CONFIG"; then
  warning 'could not prepare GHCR authentication; continuing without cache transfer'
  exit 0
fi

# Two bounded attempts. Each caller supplies its command timeout.
retry() {
  local attempt
  for attempt in 1 2; do
    if "$@" >"$log" 2>&1; then return 0; fi
    if [[ "$attempt" == 1 ]]; then sleep 2; fi
  done
  return 1
}
login() {
  printf '%s' "$GHCR_TOKEN" | timeout --kill-after=5s 30s docker login ghcr.io \
    --username "$GITHUB_ACTOR" --password-stdin
}
if [[ -n "$GHCR_TOKEN" ]] && ! retry login; then
  warning 'could not log in to GHCR; continuing without cache transfer'
  exit 0
fi

if [[ "$CACHE_MODE" == restore ]]; then
  if ! retry timeout --kill-after=5s 60s docker pull "$ref"; then
    if grep -Eqi 'manifest unknown|manifest.*not found|no such manifest' "$log"; then
      echo "no trunk at $ref: cold build"
    else
      warning "could not pull ccache trunk $ref; continuing with the local cache"
    fi
    exit 0
  fi
  # Stage beside the destination so promotion uses renames on the same filesystem.
  if ! mkdir -p -- "$(dirname -- "$path")" || ! stage=$(mktemp -d "$path.ghcr-XXXXXX"); then
    warning 'could not stage ccache restore; keeping the local cache'
    exit 0
  fi
  container="erf-ccache-$$-$RANDOM"
  if ! timeout --kill-after=5s 15s docker create --name "$container" "$ref" /trunk >"$log" 2>&1 ||
     ! timeout --kill-after=5s 60s docker cp "$container:/ccache/." "$stage/" >"$log" 2>&1; then
    warning 'could not extract ccache trunk; keeping the local cache'
    exit 0
  fi
  if [[ -e "$path" ]]; then
    backup="$stage.previous"
    if ! mv -- "$path" "$backup"; then
      warning 'could not prepare ccache replacement; keeping the local cache'
      exit 0
    fi
  fi
  if ! mv -- "$stage" "$path"; then
    warning 'could not install ccache trunk; restoring the local cache'
    exit 0
  fi
  stage=''
  echo "restored $ref: $(du -sh "$path" | cut -f1)"
  echo 'hit=true' >> "$GITHUB_OUTPUT"
  exit 0
fi

printf 'FROM scratch\nCOPY . /ccache\n' > "$scratch/Dockerfile"
if ! timeout --kill-after=5s 60s docker build --quiet -t "$ref" \
    --label "org.opencontainers.image.source=$CACHE_SOURCE" \
    --label "org.opencontainers.image.revision=$GITHUB_SHA" \
    -f "$scratch/Dockerfile" "$path" >"$log" 2>&1; then
  warning 'could not package ccache trunk; continuing without saving'
  exit 0
fi
if ! retry timeout --kill-after=5s 60s docker push --quiet "$ref"; then
  warning 'could not push ccache trunk; continuing without saving'
  exit 0
fi
echo "saved $ref: $(du -sh "$path" | cut -f1)"
# Bound the whole prune, including pagination and any retries. A later save retries it.
if ! GH_TOKEN="$GHCR_TOKEN" OWNER="$GITHUB_REPOSITORY_OWNER" OWNER_TYPE="$CACHE_OWNER_TYPE" \
    timeout --kill-after=5s 60s bash "$(dirname -- "${BASH_SOURCE[0]}")/prune.sh" erf-ccache:0 \
    >"$log" 2>&1; then
  warning 'could not prune old GHCR ccache versions; a later save will retry'
else
  cat "$log"
fi
