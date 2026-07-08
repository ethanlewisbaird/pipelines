# BAIRD Pipelines

Reusable bioinformatics pipelines for BAIRD.

## Structure

```
pipelines/
├── common/          # Shared utilities (QC, plotting, etc.)
├── scrnaseq/        # Single-cell RNA-seq pipelines
├── rnaseq/          # Bulk RNA-seq pipelines
├── spatial/         # Spatial omics pipelines
└── atac/            # ATAC-seq pipelines
```

## Usage

Pipelines are automatically cloned/cached by BAIRD workers. To use a pipeline in a job:

```bash
# In job context:
{
  "pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git",
  "script_path": "scrnaseq/pseudotime.py"
}
```

Or via CLI:
```bash
baird jobs submit --id my-job --project my-project \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "scrnaseq/pseudotime.py"}' \
  --command "python /data/ebaird/pipelines/scrnaseq/pseudotime.py --input data.h5ad"
```

## Adding a new pipeline

1. Create a directory: `pipelines/my_analysis/`
2. Add your scripts
3. Test locally
4. Push to GitHub
5. Update project conventions with `pipeline_repo`

## Versioning

- Base pipelines are in `main` branch
- Tags for stable versions: `v1.0`, `v2.0`
- BAIRD tracks commit hash for provenance
