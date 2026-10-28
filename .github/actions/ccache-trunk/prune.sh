#!/usr/bin/env bash
# Prunes superseded ccache images. Called with erf-ccache:0 after a successful push.
#
# Keeps every version that carries a tag other than sha-<commit> (development,
# latest, a release number) and the newest KEEP sha-only versions of each
# package; deletes the rest, untagged versions (superseded manifests) included.
# A failed request is reported with its response and retried; a version that
# still cannot be deleted is left for the next run, and the exit status says
# something failed.
#
# Usage: prune.sh <package>[:<keep>]...
#   erf-ccache:0   keep named variant tags; remove superseded versions
# Environment: OWNER (the repository owner), OWNER_TYPE (Organization or User),
# GH_TOKEN (packages: write on the publishing repository), KEEP (the default
# per-package keep, 10).
set -uo pipefail

: "${OWNER:?}" "${OWNER_TYPE:?}" "${GH_TOKEN:?}"
KEEP=${KEEP:-10}
readonly VERSION_ALREADY_DELETED=4

case "$OWNER_TYPE" in
  Organization) owner_api_path="orgs/$OWNER" ;;
  User) owner_api_path="users/$OWNER" ;;
  *)
    echo "unexpected OWNER_TYPE: $OWNER_TYPE" >&2
    exit 1
    ;;
esac

# Try each GitHub API request up to three times. Only successful response data goes
# to stdout; diagnostics go to stderr so callers can capture version listings.
# A concurrent pruner may have already deleted a version: report that separately.
github_api() {
  local response attempt argument
  local is_delete_request=0

  for argument in "$@"; do
    if [[ "$argument" == DELETE ]]; then
      is_delete_request=1
    fi
  done

  for attempt in 1 2 3; do
    if response=$(timeout --kill-after=5s 15s gh api "$@" 2>&1); then
      if [[ -n "$response" ]]; then
        printf '%s\n' "$response"
      fi
      return 0
    fi
    if [[ "$is_delete_request" -eq 1 && "$response" == *"HTTP 404"* ]]; then
      echo "gh api $*: version already absent (HTTP 404)" >&2
      return "$VERSION_ALREADY_DELETED"
    fi
    echo "gh api $* (attempt $attempt): $response" >&2
    if [[ "$attempt" -lt 3 ]]; then
      sleep 5
    fi
  done
  return 1
}

# Read tab-separated ID, creation time, and comma-separated tags from stdin.
# Preserve every named tag and the requested number of newest sha-only versions.
# Print IDs of the remaining sha-only versions and all untagged versions.
select_versions_to_delete() {
  local keep_sha_versions=$1

  sort -t $'\t' -k2,2r | awk -F '\t' -v keep="$keep_sha_versions" '
      NF == 0 { next }
      {
        version_id = $1
        tag_count = split($3, tags, ",")
        has_named_tag = 0
        has_sha_tag = 0

        for (i = 1; i <= tag_count; i++) {
          if (tags[i] == "") continue
          if (tags[i] ~ /^sha-/) {
            has_sha_tag = 1
          } else {
            has_named_tag = 1
          }
        }

        if (has_named_tag) next
        if (has_sha_tag) {
          sha_versions_seen++
          if (sha_versions_seen <= keep) next
        }
        print version_id
      }'
}

prune_package() {
  local package_spec=$1
  local package_name keep_sha_versions versions_api_path versions versions_to_delete
  local version_id delete_status
  local package_failed=0

  package_name=${package_spec%%:*}
  keep_sha_versions=${package_spec#*:}
  if [[ "$keep_sha_versions" == "$package_spec" ]]; then
    keep_sha_versions=$KEEP
  fi
  versions_api_path="$owner_api_path/packages/container/$package_name/versions"

  if ! versions=$(github_api --paginate "$versions_api_path?per_page=100" \
        --jq '.[] | [.id, .created_at, (.metadata.container.tags | join(","))] | @tsv'); then
    echo "$package_name: could not list versions" >&2
    return 1
  fi
  versions_to_delete=$(printf '%s\n' "$versions" | select_versions_to_delete "$keep_sha_versions")

  while read -r version_id; do
    if [[ -z "$version_id" ]]; then
      continue
    fi
    if github_api --method DELETE --silent "$versions_api_path/$version_id"; then
      echo "$package_name: deleted version $version_id"
    else
      delete_status=$?
      if [[ "$delete_status" -eq "$VERSION_ALREADY_DELETED" ]]; then
        echo "$package_name: version $version_id was already gone"
      else
        echo "$package_name: could not delete version $version_id" >&2
        package_failed=1
      fi
    fi
  done <<< "$versions_to_delete"

  return "$package_failed"
}

had_failures=0
for package_spec in "$@"; do
  if ! prune_package "$package_spec"; then
    had_failures=1
  fi
done
exit "$had_failures"
