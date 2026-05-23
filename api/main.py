import json
import os
import subprocess
import tempfile
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import JSONResponse

app = FastAPI(
    title="Gene Annotation HMM API",
    description=(
        "Annotates genome sequences using pre-trained HMM models. "
        "POST a FASTA file and a genome name; the server decodes it with Viterbi "
        "and returns per-chromosome gene structures (exons, introns, boundaries)."
    ),
    version="1.0.0",
)

MODELS_DIR  = Path(os.getenv("MODELS_DIR",  "/app/models"))
ANNOTATE_BIN = os.getenv("ANNOTATE_BIN", "/usr/local/bin/annotate")


def _available_genomes() -> list[str]:
    if not MODELS_DIR.exists():
        return []
    return sorted(
        p.stem.removesuffix("_model")
        for p in MODELS_DIR.glob("*_model.json")
    )


@app.get("/health", summary="Health check")
def health():
    return {
        "status": "ok",
        "annotate_bin": ANNOTATE_BIN,
        "models_dir": str(MODELS_DIR),
        "available_genomes": _available_genomes(),
    }


@app.get("/genomes", summary="List available genome models")
def list_genomes():
    """Returns the names of all pre-trained genome models stored in the container."""
    return {"genomes": _available_genomes()}


@app.post("/annotate", summary="Annotate a genome FASTA file")
async def annotate(
    genome: str = Form(
        ...,
        description="Genome profile name matching a pre-trained model (e.g. 'yeast').",
    ),
    fna: UploadFile = File(
        ...,
        description="FASTA genome file (.fna / .fa / .fasta). "
                    "All chromosomes present will be decoded.",
    ),
):
    """
    Runs Viterbi decoding on the uploaded FASTA using the named genome model.

    Returns a JSON object with:
    - `genome`: model name used
    - `chromosomes`: list of decoded chromosomes, each with `genes` (exons + introns)
    - `summary`: total gene and chromosome counts
    """
    model_path = MODELS_DIR / f"{genome}_model.json"
    if not model_path.exists():
        available = _available_genomes()
        raise HTTPException(
            status_code=404,
            detail=f"No trained model for '{genome}'. Available: {available}",
        )

    content = await fna.read()
    with tempfile.NamedTemporaryFile(suffix=".fna", delete=False) as tmp:
        tmp.write(content)
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            [ANNOTATE_BIN, "--model", str(model_path), "--fasta", tmp_path],
            capture_output=True,
            text=True,
            timeout=600,
        )
    except subprocess.TimeoutExpired:
        raise HTTPException(status_code=504, detail="Annotation timed out (>10 min)")
    finally:
        os.unlink(tmp_path)

    if result.returncode != 0:
        raise HTTPException(
            status_code=500,
            detail=f"Annotation binary failed: {result.stderr.strip()}",
        )

    try:
        annotation = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Could not parse annotator output: {exc}",
        )

    return JSONResponse(content=annotation)
