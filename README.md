# RAPIDS → AlphaFold2/AlphaFold DB → GROMACS Workflow

This project runs a small workflow for:

```text
single-cell analysis [RAPIDS-singlecell]
→ marker gene selection
→ UniProt protein FASTA retrieval
→ AlphaFold DB / NVIDIA AlphaFold2 API
→ GROMACS energy minimization
```

CUDA-Q is not part of this version.

## Required input data

The default input is:

```text
data/pbmc3k.h5ad
```

This is an AnnData h5ad file containing a single-cell gene expression matrix. You can download the example dataset with:

```bash
python src/download_singlecell.py --config config.yaml
```

You may replace it with another `.h5ad` dataset by editing:

```yaml
single_cell:
  input_h5ad: data/your_dataset.h5ad
```

## Recommended WSL2 setup

### Option A: Python venv

```bash
cd ~/projects/nim_blueprint

conda deactivate  # repeat until only the normal shell remains if needed

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements.txt
```

Install RAPIDS-singlecell separately according to your CUDA version. For CUDA 12:

```bash
python -m pip install "rapids-singlecell-cu12[rapids]" --extra-index-url=https://pypi.nvidia.com
```

Test GPU visibility:

```bash
python -c "import cupy as cp; print(cp.cuda.runtime.getDeviceCount())"
```

If CuPy reports `libnvrtc.so.12` missing, install NVRTC in the same environment or system:

```bash
conda install -c conda-forge cuda-nvrtc
```

### Option B: Conda environment

```bash
conda create -n nim-blueprint python=3.11 -y
conda activate nim-blueprint
pip install -r requirements.txt
pip install "rapids-singlecell-cu12[rapids]" --extra-index-url=https://pypi.nvidia.com
```

## GROMACS setup

GROMACS can be used either as a local command or through a Docker wrapper.

### Local / conda GROMACS

Install with conda-forge:

```bash
conda create -n gromacs -c conda-forge gromacs -y
conda activate gromacs
gmx --version
```

If you want the Python workflow to call this GROMACS from another shell, set the full path in `config.yaml`:

```yaml
gromacs:
  gmx_bin: /home/user/miniconda3/envs/gromacs/bin/gmx
```

### Docker GROMACS

Make the wrapper executable:

```bash
chmod +x scripts/gmx_docker.sh
```

Set `config.yaml`:

```yaml
gromacs:
  gmx_bin: scripts/gmx_docker.sh
```

Optionally choose an image:

```bash
export GMX_IMAGE="gromacs/gromacs:latest"
```

Then run the workflow normally. The wrapper mounts the GROMACS working directory into the container.

## AlphaFold2 settings

The workflow checks AlphaFold DB first by UniProt accession. If no structure is available, it can call NVIDIA AlphaFold2 API.

For NVIDIA API usage:

```bash
export NVIDIA_API_KEY="nvapi-..."
```

The endpoint and polling options are configured in `config.yaml`.

If the hosted API returns HTTP 500/504, use the AlphaFold DB-first mode for known proteins. This avoids unnecessary API calls when a public structure already exists.

## Run workflow step by step

```bash
source .venv/bin/activate

python src/download_singlecell.py --config config.yaml
python src/rapids_singlecell_.py --config config.yaml
python src/gene_to_fasta.py --config config.yaml
python src/alphafold2_api.py --config config.yaml
python src/gromacs.py --config config.yaml
python src/report.py --config config.yaml
```

## Run full pipeline

```bash
python pipeline.py
```

## Main outputs

```text
outputs/single_cell/marker_genes.csv
outputs/single_cell/top_marker_genes.csv
outputs/fasta/*.fasta
outputs/folding/target.pdb
outputs/gromacs/em.gro
outputs/gromacs/em.edr
outputs/gromacs/em.log
outputs/report.md
```

## Notes

This version stops at protein structure preparation and GROMACS energy minimization. For drug discovery, the next stages would typically be docking, ligand screening, longer MD, and binding energy analysis.
