import cors from "cors";
import express from "express";
import { execFile } from "node:child_process";
import { randomUUID } from "node:crypto";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import multer from "multer";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const repoRoot = path.resolve(__dirname, "..");
const localDataDir = path.join(__dirname, "local_data");
const uploadDir = path.join(localDataDir, "uploads");
const binDir = path.join(localDataDir, "bin");
const predictorBin = path.join(binDir, "hmm_predict_fna");

fs.mkdirSync(uploadDir, { recursive: true });
fs.mkdirSync(binDir, { recursive: true });

const app = express();
const upload = multer({ dest: uploadDir });
const runs = new Map();

app.use(cors());
app.use(express.json());

function execFilePromise(command, args, options = {}) {
  return new Promise((resolve, reject) => {
    execFile(command, args, { maxBuffer: 1024 * 1024 * 256, ...options }, (error, stdout, stderr) => {
      if (error) {
        error.stderr = stderr;
        reject(error);
        return;
      }
      resolve({ stdout, stderr });
    });
  });
}

async function ensurePredictorBuilt() {
  const sources = [
    "src/tools/predict_fna.cpp",
    "src/decoding/Viterbi.cpp",
    "src/decoding/Forward_Backward.cpp",
    "src/model/Transition_Model.cpp",
    "src/genome_profiles/Genome_Profile.cpp",
    "src/parsers/FNA_Parser.cpp",
    "src/parsers/GFF_Parser.cpp",
    "src/model/Emission_Model.cpp",
  ].map((source) => path.join(repoRoot, source));

  if (fs.existsSync(predictorBin)) {
    const binaryMtime = fs.statSync(predictorBin).mtimeMs;
    const newestSourceMtime = Math.max(...sources.map((source) => fs.statSync(source).mtimeMs));
    if (binaryMtime >= newestSourceMtime) return;
  }

  await execFilePromise(
    "clang++",
    [
      "-std=c++17",
      "-Isrc",
      "-I/opt/homebrew/include",
      "src/tools/predict_fna.cpp",
      "src/decoding/Viterbi.cpp",
      "src/decoding/Forward_Backward.cpp",
      "src/model/Transition_Model.cpp",
      "src/genome_profiles/Genome_Profile.cpp",
      "src/parsers/FNA_Parser.cpp",
      "src/parsers/GFF_Parser.cpp",
      "src/model/Emission_Model.cpp",
      "-o",
      predictorBin,
    ],
    { cwd: repoRoot }
  );
}

function makeRunName(fileName) {
  const stem = path.basename(fileName).replace(/\.(fna|fa|fasta)$/i, "");
  const stamp = new Date().toISOString().replace(/[-:]/g, "").slice(0, 15);
  return `${stem}_${stamp}`;
}

function makeRunRecord(file, predictionResult, elapsedMs) {
  const id = randomUUID();
  const run = {
    id,
    name: makeRunName(file.originalname),
    fileName: file.originalname,
    date: new Date().toLocaleDateString("en-US", {
      month: "short",
      day: "numeric",
      year: "numeric",
    }),
    status: "done",
    elapsedMs,
    summary: predictionResult.summary,
    scaffolds: predictionResult.scaffolds,
    predictions: predictionResult.predictions,
    confidenceByScaffold: predictionResult.confidenceByScaffold,
  };
  runs.set(id, run);
  return run;
}

function escapeAttribute(value) {
  return String(value).replace(/[;\t\n\r]/g, "_");
}

function runToGff3(run) {
  const lines = ["##gff-version 3"];
  for (const prediction of run.predictions) {
    const strand = prediction.strand === "-" ? "-" : "+";
    lines.push(
      [
        prediction.scaffold,
        "HMMGenePredictor",
        "gene",
        prediction.start,
        prediction.end,
        ".",
        strand,
        ".",
        `ID=${escapeAttribute(prediction.id)}`,
      ].join("\t")
    );
    prediction.exons.forEach((exon, index) => {
      lines.push(
        [
          prediction.scaffold,
          "HMMGenePredictor",
          "CDS",
          exon.start,
          exon.end,
          ".",
          strand,
          ".",
          `ID=${escapeAttribute(prediction.id)}.cds${index + 1};Parent=${escapeAttribute(prediction.id)}`,
        ].join("\t")
      );
    });
  }
  return `${lines.join("\n")}\n`;
}

function runToCsv(run) {
  const rows = ["id,scaffold,strand,start,end,length,exon_count,intron_count"];
  for (const prediction of run.predictions) {
    rows.push(
      [
        prediction.id,
        prediction.scaffold,
        prediction.strand ?? "+",
        prediction.start,
        prediction.end,
        prediction.end - prediction.start + 1,
        prediction.exons.length,
        prediction.introns.length,
      ].join(",")
    );
  }
  return `${rows.join("\n")}\n`;
}

function runToBed(run) {
  const lines = [];
  for (const prediction of run.predictions) {
    lines.push(
      [
        prediction.scaffold,
        prediction.start - 1,
        prediction.end,
        prediction.id,
        0,
        prediction.strand ?? "+",
      ].join("\t")
    );
  }
  return `${lines.join("\n")}\n`;
}

function runToFasta(run) {
  const lines = [];
  for (const prediction of run.predictions) {
    lines.push(`>${prediction.id} ${prediction.scaffold}:${prediction.start}-${prediction.end}`);
    lines.push(prediction.sequence || "");
  }
  return `${lines.join("\n")}\n`;
}

app.get("/api/health", (_req, res) => {
  res.json({ ok: true, mode: "local" });
});

app.post("/api/predict", upload.single("file"), async (req, res) => {
  if (!req.file) {
    res.status(400).json({ error: "FASTA file is required." });
    return;
  }

  const started = Date.now();
  try {
    await ensurePredictorBuilt();
    const { stdout } = await execFilePromise(
      predictorBin,
      [
        "--fna",
        req.file.path,
        "--profile",
        path.join(repoRoot, "src/genome_profiles/yeast.json"),
      ],
      { cwd: repoRoot }
    );
    const predictionResult = JSON.parse(stdout);
    const run = makeRunRecord(req.file, predictionResult, Date.now() - started);
    res.json(run);
  } catch (error) {
    res.status(500).json({
      error: "Prediction failed.",
      detail: error.stderr || error.message,
    });
  }
});

app.get("/api/export/:runId/:format", (req, res) => {
  const run = runs.get(req.params.runId);
  if (!run) {
    res.status(404).json({ error: "Run not found." });
    return;
  }

  const format = req.params.format;
  const exporters = {
    bed: ["text/plain", runToBed],
    csv: ["text/csv", runToCsv],
    fasta: ["text/plain", runToFasta],
    gff3: ["text/plain", runToGff3],
  };
  const exporter = exporters[format];
  if (!exporter) {
    res.status(400).json({ error: "Unsupported export format." });
    return;
  }

  const [contentType, render] = exporter;
  res.setHeader("Content-Type", contentType);
  res.setHeader("Content-Disposition", `attachment; filename="${run.name}.${format}"`);
  res.send(render(run));
});

const port = Number(process.env.PORT || 5174);
const httpServer = app.listen(port, () => {
  console.log(`Local HMM API listening on http://localhost:${port}`);
});

httpServer.on("error", (error) => {
  if (error.code === "EADDRINUSE") {
    console.error(`Port ${port} is already in use. Stop the other process or set PORT.`);
  } else {
    console.error(error);
  }
  process.exit(1);
});
