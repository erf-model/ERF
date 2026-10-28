#!/usr/bin/env bash
# Linux only: Docker and GNU coreutils are available on the Ubuntu CI runners.
# Invalid inputs fail; cache transport failures warn and leave CI free to continue.
set -euo pipefail

error() {
  echo "::error::$*" >&2
  exit 2
}

warning() {
  if [[ -f ${transfer_log:-} ]]; then
    tail -n 20 "$transfer_log" >&2
  fi
  echo "::warning::$*" >&2
}

case "$CACHE_MODE" in
  restore|save) ;;
  *) error 'mode must be restore or save' ;;
esac
if [[ "$RUNNER_OS" != Linux ]]; then
  error 'ccache-trunk requires a Linux runner'
fi
if [[ "$CACHE_MODE" == save && -z "$GHCR_TOKEN" ]]; then
  error 'token is required when saving'
fi

owner=$(printf '%s' "$GITHUB_REPOSITORY_OWNER" | tr '[:upper:]' '[:lower:]')
tag=$(printf '%s' "$CACHE_KEY" | sed -E 's/[^A-Za-z0-9_.-]+/-/g' | tr '[:upper:]' '[:lower:]')
if [[ ! "$tag" =~ ^[a-z0-9_][a-z0-9_.-]{0,127}$ ]]; then
  error 'key must produce a valid GHCR tag'
fi
image_ref="ghcr.io/$owner/erf-ccache:$tag"

case "$CACHE_PATH" in
  \~) cache_dir=$HOME ;;
  \~/*) cache_dir="$HOME/${CACHE_PATH:2}" ;;
  '') error 'path must name a dedicated ccache directory' ;;
  *) cache_dir=$CACHE_PATH ;;
esac
cache_dir=$(realpath -m -- "$cache_dir")
if [[ "$cache_dir" == / || "$cache_dir" == "$HOME" ]]; then
  error 'path must name a dedicated ccache directory'
fi
if [[ -e "$cache_dir" && ! -d "$cache_dir" ]]; then
  error 'cache path is not a directory'
fi

if [[ "$CACHE_MODE" == restore ]]; then
  echo 'hit=false' >> "$GITHUB_OUTPUT"
elif [[ ! -d "$cache_dir" || -z $(find "$cache_dir" -type f -print -quit) ]]; then
  echo "no ccache files in $cache_dir: skip save"
  exit 0
fi

if ! transfer_dir=$(mktemp -d "$RUNNER_TEMP/erf-ccache.XXXXXX"); then
  warning 'could not prepare GHCR cache transfer; continuing without it'
  exit 0
fi
transfer_log="$transfer_dir/transfer.log"
restore_container=''
staged_cache_dir=''
previous_cache_dir=''

# restore_cache replaces the cache with two renames:
#   cache_dir -> previous_cache_dir, then staged_cache_dir -> cache_dir.
# On exit, a missing cache_dir means promotion did not finish: put the old cache back.
# If cache_dir exists, promotion succeeded and the old directory can be removed.
# These paths are shared with restore_cache so cleanup also works after interruption.
cleanup() {
  if [[ -n "$restore_container" ]]; then
    timeout --kill-after=5s 10s docker rm "$restore_container" >/dev/null 2>&1 || true
  fi

  if [[ -n "$previous_cache_dir" && -e "$previous_cache_dir" ]]; then
    if [[ ! -e "$cache_dir" ]]; then
      if ! mv -- "$previous_cache_dir" "$cache_dir"; then
        echo "::warning::original cache remains at $previous_cache_dir" >&2
      fi
    else
      rm -rf -- "$previous_cache_dir" || true
    fi
  fi

  if [[ -n "$staged_cache_dir" ]]; then
    rm -rf -- "$staged_cache_dir" || true
  fi
  rm -rf -- "$transfer_dir" || true
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

# Keep credentials out of the runner's shared Docker config and remove them on exit.
export DOCKER_CONFIG="$transfer_dir/docker"
if ! mkdir -p "$DOCKER_CONFIG"; then
  warning 'could not prepare GHCR authentication; continuing without cache transfer'
  exit 0
fi

# Two bounded attempts. Each caller supplies its command timeout.
retry() {
  local attempt
  for attempt in 1 2; do
    if "$@" >"$transfer_log" 2>&1; then
      return 0
    fi
    if [[ "$attempt" == 1 ]]; then
      sleep 2
    fi
  done
  return 1
}

login_to_registry() {
  printf '%s' "$GHCR_TOKEN" | timeout --kill-after=5s 30s docker login ghcr.io \
    --username "$GITHUB_ACTOR" --password-stdin
}

restore_cache() {
  if ! retry timeout --kill-after=5s 60s docker pull "$image_ref"; then
    if grep -Eqi 'manifest unknown|manifest.*not found|no such manifest' "$transfer_log"; then
      echo "no trunk at $image_ref: cold build"
    else
      warning "could not pull ccache trunk $image_ref; continuing with the local cache"
    fi
    return 0
  fi

  # Stage beside the destination so promotion uses renames on the same filesystem.
  if ! mkdir -p -- "$(dirname -- "$cache_dir")"; then
    warning 'could not stage ccache restore; keeping the local cache'
    return 0
  fi
  if ! staged_cache_dir=$(mktemp -d "$cache_dir.ghcr-XXXXXX"); then
    warning 'could not stage ccache restore; keeping the local cache'
    return 0
  fi

  restore_container="erf-ccache-$$-$RANDOM"
  if ! timeout --kill-after=5s 15s docker create --name "$restore_container" "$image_ref" /trunk \
      >"$transfer_log" 2>&1; then
    warning 'could not extract ccache trunk; keeping the local cache'
    return 0
  fi
  if ! timeout --kill-after=5s 60s docker cp "$restore_container:/ccache/." "$staged_cache_dir/" \
      >"$transfer_log" 2>&1; then
    warning 'could not extract ccache trunk; keeping the local cache'
    return 0
  fi

  # Keep the old directory until EXIT cleanup can confirm promotion succeeded.
  if [[ -e "$cache_dir" ]]; then
    previous_cache_dir="$staged_cache_dir.previous"
    if ! mv -- "$cache_dir" "$previous_cache_dir"; then
      warning 'could not prepare ccache replacement; keeping the local cache'
      return 0
    fi
  fi
  if ! mv -- "$staged_cache_dir" "$cache_dir"; then
    warning 'could not install ccache trunk; restoring the local cache'
    return 0
  fi
  staged_cache_dir=''
  echo "restored $image_ref: $(du -sh "$cache_dir" | cut -f1)"
  echo 'hit=true' >> "$GITHUB_OUTPUT"
}

save_cache() {
  local prune_script

  printf 'FROM scratch\nCOPY . /ccache\n' > "$transfer_dir/Dockerfile"
  if ! timeout --kill-after=5s 60s docker build --quiet -t "$image_ref" \
      --label "org.opencontainers.image.source=$CACHE_SOURCE" \
      --label "org.opencontainers.image.revision=$GITHUB_SHA" \
      -f "$transfer_dir/Dockerfile" "$cache_dir" >"$transfer_log" 2>&1; then
    warning 'could not package ccache trunk; continuing without saving'
    return 0
  fi
  if ! retry timeout --kill-after=5s 60s docker push --quiet "$image_ref"; then
    warning 'could not push ccache trunk; continuing without saving'
    return 0
  fi
  echo "saved $image_ref: $(du -sh "$cache_dir" | cut -f1)"

  # Bound the whole prune, including pagination and any retries. A later save retries it.
  prune_script="$(dirname -- "${BASH_SOURCE[0]}")/prune.sh"
  if ! GH_TOKEN="$GHCR_TOKEN" OWNER="$GITHUB_REPOSITORY_OWNER" OWNER_TYPE="$CACHE_OWNER_TYPE" \
      timeout --kill-after=5s 60s bash "$prune_script" erf-ccache:0 >"$transfer_log" 2>&1; then
    warning 'could not prune old GHCR ccache versions; a later save will retry'
  else
    cat "$transfer_log"
  fi
}

if [[ -n "$GHCR_TOKEN" ]]; then
  if ! retry login_to_registry; then
    warning 'could not log in to GHCR; continuing without cache transfer'
    exit 0
  fi
fi

case "$CACHE_MODE" in
  restore) restore_cache ;;
  save) save_cache ;;
esac
