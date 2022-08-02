Environment variables that help to set are:
- `export NCBI_API_KEY="api key"` to increase rate limit with NCBI
- `export TMPDIR=/ebs/tmpfs/tmp` to have more tmp file space for Python (after running, e.g., `aws-tmpfs-vol-create 100`)
