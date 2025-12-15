# Supply chain: signing + attestations

crispr-gpu publishes a CPU Docker image to GHCR and signs it **keylessly** with **cosign** in GitHub Actions.

## Important warning

Do **not** use `:latest` for scientific claims. `:latest` is for convenience only.

For anything citeable/reproducible, pin either:
- a release tag (`:vX.Y.Z`), or
- the immutable digest (`@sha256:...`) — preferred.

If you are verifying a specific released tag, prefer the canonical digest and copy/paste verification commands in that tag’s GitHub Release. This document tracks the current recommended verification flags and may evolve without changing older release artifacts.

## Verify a release image signature (keyless)

Set these:
```bash
IMAGE_WITH_DIGEST="ghcr.io/omniscoder/crispr-gpu-cpu@sha256:..."
```

Verify the signature for the digest:
```bash
cosign verify \
  --certificate-identity-regexp '^https://github.com/omniscoder/crispr-gpu/.github/workflows/docker-publish\.yml@refs/(tags/v[0-9]+\.[0-9]+\.[0-9]+|heads/master)$' \
  --certificate-oidc-issuer 'https://token.actions.githubusercontent.com' \
  "${IMAGE_WITH_DIGEST}"
```

If you want *strict release-only* verification, tighten the identity to `@refs/tags/vX.Y.Z` instead of allowing `@refs/heads/master`.

## Verify an SBOM attestation exists (same digest)

```bash
cosign verify-attestation \
  --type spdxjson \
  --certificate-identity-regexp '^https://github.com/omniscoder/crispr-gpu/.github/workflows/docker-publish\.yml@refs/(tags/v[0-9]+\.[0-9]+\.[0-9]+|heads/master)$' \
  --certificate-oidc-issuer 'https://token.actions.githubusercontent.com' \
  "${IMAGE_WITH_DIGEST}"
```

## Embed image provenance into `report.json`

If you run the demo inside the container and want the report to be self-describing, inject the digest (and optional cosign signature reference) into the container environment:
```bash
IMAGE="ghcr.io/omniscoder/crispr-gpu-cpu"
DIGEST="sha256:..."   # immutable digest, not a tag
SIG="$(cosign triangulate --type signature "${IMAGE}@${DIGEST}")"
docker run --rm \
  -e CRISPR_GPU_IMAGE_REPOSITORY="${IMAGE}" \
  -e CRISPR_GPU_IMAGE_DIGEST="${DIGEST}" \
  -e CRISPR_GPU_COSIGN_SIGNATURE="${SIG}" \
  -v "$(pwd)/reports_docker:/out" \
  "${IMAGE}@${DIGEST}"
```
