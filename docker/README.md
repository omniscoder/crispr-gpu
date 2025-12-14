# Docker

These images are intended for artifact evaluation: run a deterministic demo + synthetic benchmark and emit a versioned report bundle.

## CPU image

Build:
```bash
docker build -f docker/Dockerfile.cpu -t crispr-gpu:cpu .
```

Run (writes into `./reports_docker` on the host):
```bash
mkdir -p reports_docker
docker run --rm -v "$(pwd)/reports_docker:/out" crispr-gpu:cpu
```

## CUDA image

Build:
```bash
docker build -f docker/Dockerfile.cuda -t crispr-gpu:cuda .
```

Run (requires NVIDIA Container Toolkit):
```bash
mkdir -p reports_docker
docker run --rm --gpus all -v "$(pwd)/reports_docker:/out" crispr-gpu:cuda
```

Outputs:
- `/out/report.json` (schema: `schemas/report.v1.json`)
- `/out/report.md`

