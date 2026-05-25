import {
  BarChart3,
  ChevronDown,
  ChevronRight,
  Copy,
  Download,
  FileArchive,
  FileCode2,
  FileText,
  FlaskConical,
  Folder,
  FolderOpen,
  MoreHorizontal,
  Play,
  Plus,
  Search,
  UploadCloud,
} from "lucide-react";
import { useMemo, useRef, useState } from "react";

const savedRuns = [
  {
    id: "yeast_test_01",
    name: "yeast_test_01",
    fileName: "yeast_chr01.fna",
    date: "May 24, 2026",
    status: "done",
    totalBases: 230218,
    scaffolds: 1,
    predictedGenes: 112,
    predictedExons: 224,
    predictedIntrons: 112,
    scaffold: "chr01",
  },
  {
    id: "chr1_pre_alpha",
    name: "chr1_pre_alpha",
    fileName: "human_chr1_sample.fna",
    date: "May 20, 2026",
    status: "done",
    totalBases: 185420,
    scaffolds: 1,
    predictedGenes: 84,
    predictedExons: 173,
    predictedIntrons: 89,
    scaffold: "chr1",
  },
  {
    id: "mouse_chrX_run",
    name: "mouse_chrX_run",
    fileName: "mouse_chrX_region.fna",
    date: "May 15, 2026",
    status: "failed",
    totalBases: 0,
    scaffolds: 0,
    predictedGenes: 0,
    predictedExons: 0,
    predictedIntrons: 0,
    scaffold: "chrX",
  },
  {
    id: "drosophila_04",
    name: "drosophila_04",
    fileName: "dmel_scaffold_04.fna",
    date: "May 12, 2026",
    status: "done",
    totalBases: 298771,
    scaffolds: 3,
    predictedGenes: 139,
    predictedExons: 281,
    predictedIntrons: 142,
    scaffold: "scaffold_04",
  },
];

const predictionsByRun = {
  yeast_test_01: [
    {
      id: "pred_0001",
      scaffold: "chr01",
      start: 14600,
      end: 16400,
      exons: [
        [14600, 15450],
        [15820, 16400],
      ],
      introns: [[15451, 15819]],
      sequence:
        "ATGGCGTCACTACGGGCATAGCTACGTAGCTACGTACTAGCTACGCTGACTGACTGATCGATCGATCGTAGTCGTAGCTGATCGTAGCTGACTGATCGATCGTAGCTGACTGATCGTAGCTAATAG",
    },
    {
      id: "pred_0002",
      scaffold: "chr01",
      start: 17600,
      end: 19100,
      exons: [[17600, 19100]],
      introns: [],
      sequence:
        "ATGCTAGCTACGATCGATCGATGCTAGCTAGCTACGATCGTACGTACGATCGATCGTAGCTAGCTACGATCGTAGCTAGCTAGCTGACTAG",
    },
    {
      id: "pred_0003",
      scaffold: "chr01",
      start: 21200,
      end: 24500,
      exons: [
        [21200, 21930],
        [22480, 23110],
        [23550, 24120],
        [24280, 24500],
      ],
      introns: [
        [21931, 22479],
        [23111, 23549],
        [24121, 24279],
      ],
      sequence:
        "ATGACGTACTGACTGATCGTACGTACTGATCGATCGTAGCTGACTAGCTAGCTAGCTGACTGATCGTAGCTAGCTAGCTGATCGATCGTAA",
    },
  ],
  chr1_pre_alpha: [
    {
      id: "pred_0001",
      scaffold: "chr1",
      start: 8020,
      end: 10990,
      exons: [
        [8020, 8900],
        [9480, 10990],
      ],
      introns: [[8901, 9479]],
      sequence:
        "ATGCGTACGATCGATCGTACGATCGTAGCTAGCTAGCTACGATCGTAGCTAGCATCGATCGATCGTAGCTAGCTACGATCGATCGTAA",
    },
    {
      id: "pred_0002",
      scaffold: "chr1",
      start: 44010,
      end: 46240,
      exons: [[44010, 46240]],
      introns: [],
      sequence:
        "ATGATCGATCGTAGCTAGCTAGCTGACTGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTTGA",
    },
  ],
  drosophila_04: [
    {
      id: "pred_0001",
      scaffold: "scaffold_04",
      start: 1200,
      end: 3160,
      exons: [
        [1200, 1840],
        [2110, 3160],
      ],
      introns: [[1841, 2109]],
      sequence:
        "ATGCTGACTGACTGATCGATCGTAGCTAGCTACGATCGATCGTAGCTAGCATCGTAGCTGACTGATCGATCGTAGCTAGCTAA",
    },
  ],
};

const emptyUpload = {
  name: "No file selected",
  size: "--",
  scaffolds: "--",
  totalBases: "--",
};

function App() {
  const [activePage, setActivePage] = useState("run");
  const [runsExpanded, setRunsExpanded] = useState(true);
  const [selectedRunId, setSelectedRunId] = useState("yeast_test_01");
  const [selectedPredictionId, setSelectedPredictionId] = useState("pred_0001");
  const [uploadedFile, setUploadedFile] = useState({
    name: "yeast_chr01.fna",
    size: "1.2 MB",
    scaffolds: 1,
    totalBases: 230218,
  });
  const [runStatus, setRunStatus] = useState("Ready");
  const [elapsedTime, setElapsedTime] = useState("0s");
  const fileInputRef = useRef(null);

  const selectedRun = useMemo(
    () => savedRuns.find((run) => run.id === selectedRunId) ?? savedRuns[0],
    [selectedRunId]
  );

  const predictions = predictionsByRun[selectedRun.id] ?? [];

  const selectedPrediction =
    predictions.find((prediction) => prediction.id === selectedPredictionId) ??
    predictions[0];

  function handleNewAnalysis() {
    setActivePage("run");
    setSelectedRunId("");
    setSelectedPredictionId("");
    setUploadedFile(emptyUpload);
    setRunStatus("Not Started");
    setElapsedTime("0s");
  }

  function openSavedRun(runId) {
    const runPredictions = predictionsByRun[runId] ?? [];
    setSelectedRunId(runId);
    setSelectedPredictionId(runPredictions[0]?.id ?? "");
    setActivePage("results");
  }

  function handleRunPrediction() {
    setRunStatus("Running");
    setElapsedTime("1s");
    window.setTimeout(() => {
      setRunStatus("Complete");
      setElapsedTime("3s");
      setSelectedRunId("yeast_test_01");
      setSelectedPredictionId("pred_0001");
      setActivePage("results");
    }, 900);
  }

  async function handleFileChange(event) {
    const file = event.target.files?.[0];
    if (!file) return;

    const text = await file.text();
    const sequenceCount = (text.match(/^>/gm) ?? []).length || 1;
    const totalBases = text
      .split("\n")
      .filter((line) => !line.startsWith(">"))
      .join("")
      .replace(/[^ACGTNacgtn]/g, "").length;

    setUploadedFile({
      name: file.name,
      size: formatBytes(file.size),
      scaffolds: sequenceCount,
      totalBases,
    });
    setRunStatus("Ready");
    setElapsedTime("0s");
  }

  return (
    <div className="app-shell">
      <aside className="sidebar">
        <div className="brand">
          <div className="brand-mark">H</div>
          <div>
            <h1>HMM Gene Predictor</h1>
            <p>v2.4.0-alpha</p>
          </div>
        </div>

        <button className="new-analysis-button" onClick={handleNewAnalysis}>
          <Plus size={17} />
          New Analysis
        </button>

        <nav className="nav-list" aria-label="Main navigation">
          <button
            className={activePage === "run" ? "nav-item active" : "nav-item"}
            onClick={() => setActivePage("run")}
          >
            <FlaskConical size={18} />
            Run Prediction
          </button>
          <button
            className={
              activePage === "results" ? "nav-item active" : "nav-item"
            }
            onClick={() => setActivePage("results")}
          >
            <BarChart3 size={18} />
            Prediction Results
          </button>
        </nav>

        <div className="saved-runs">
          <button
            className="saved-runs-header"
            onClick={() => setRunsExpanded((expanded) => !expanded)}
          >
            {runsExpanded ? <ChevronDown size={16} /> : <ChevronRight size={16} />}
            {runsExpanded ? <FolderOpen size={17} /> : <Folder size={17} />}
            <span>Saved Runs</span>
          </button>

          {runsExpanded && (
            <div className="saved-run-list">
              {savedRuns.map((run) => (
                <button
                  className={
                    selectedRunId === run.id
                      ? "saved-run active"
                      : "saved-run"
                  }
                  key={run.id}
                  onClick={() => openSavedRun(run.id)}
                >
                  <FileText size={14} />
                  <span className="saved-run-copy">
                    <span>{run.name}</span>
                    <small>{run.date}</small>
                  </span>
                  <span className={`status-pill ${run.status}`}>
                    {run.status === "done" ? "Done" : "Fail"}
                  </span>
                </button>
              ))}
            </div>
          )}
        </div>

        <div className="sidebar-footer">
          <button>
            <FileText size={17} />
            Documentation
          </button>
          <button>
            <MoreHorizontal size={17} />
            Support
          </button>
        </div>
      </aside>

      <main className="content">
        {activePage === "run" ? (
          <RunPredictionPage
            elapsedTime={elapsedTime}
            fileInputRef={fileInputRef}
            onFileChange={handleFileChange}
            onRunPrediction={handleRunPrediction}
            runStatus={runStatus}
            uploadedFile={uploadedFile}
          />
        ) : (
          <PredictionResultsPage
            predictions={predictions}
            selectedPrediction={selectedPrediction}
            selectedRun={selectedRun}
            selectedPredictionId={selectedPredictionId}
            setSelectedPredictionId={setSelectedPredictionId}
          />
        )}
      </main>
    </div>
  );
}

function RunPredictionPage({
  elapsedTime,
  fileInputRef,
  onFileChange,
  onRunPrediction,
  runStatus,
  uploadedFile,
}) {
  return (
    <section className="page">
      <PageHeader
        title="Run Prediction"
        subtitle="Ab-initio eukaryotic gene prediction from FASTA sequence"
      />

      <div className="card upload-card">
        <SectionTitle icon={<UploadCloud size={20} />} title="Input Genome" />
        <div
          className="dropzone"
          onClick={() => fileInputRef.current?.click()}
          role="button"
          tabIndex={0}
        >
          <UploadCloud size={34} />
          <p>Drag and drop FASTA file here or click to browse.</p>
          <button type="button">Choose FNA File</button>
          <small>Supported formats: .fna, .fa, .fasta</small>
          <input
            accept=".fna,.fa,.fasta"
            onChange={onFileChange}
            ref={fileInputRef}
            type="file"
          />
        </div>

        <div className="file-summary">
          <div className="file-row">
            <div>
              <FileText size={18} />
              <strong>{uploadedFile.name}</strong>
            </div>
            <span>{uploadedFile.size}</span>
          </div>
          <div className="summary-grid two">
            <Metric label="Sequences/Scaffolds" value={uploadedFile.scaffolds} />
            <Metric label="Total Bases" value={formatNumber(uploadedFile.totalBases)} />
          </div>
        </div>
      </div>

      <div className="card run-card">
        <SectionTitle icon={<Play size={20} />} title="Prediction Run" />
        <div className="run-status">
          <div>
            <p>
              <span className={`status-dot ${runStatus.toLowerCase()}`} />
              <strong>Status:</strong> {runStatus}
            </p>
            <p>
              <strong>Input:</strong> {uploadedFile.name}
            </p>
            {runStatus === "Running" && (
              <p className="running-copy">Decoding sequence with HMM...</p>
            )}
          </div>
          <div className="elapsed">
            <span>Elapsed Time</span>
            <strong>{elapsedTime}</strong>
          </div>
        </div>

        <button className="primary-action" onClick={onRunPrediction}>
          <Play size={17} />
          Run Prediction
        </button>
      </div>
    </section>
  );
}

function PredictionResultsPage({
  predictions,
  selectedPrediction,
  selectedRun,
  selectedPredictionId,
  setSelectedPredictionId,
}) {
  const [scaffold, setScaffold] = useState(selectedRun.scaffold);
  const [start, setStart] = useState(14600);
  const [end, setEnd] = useState(14699);
  const [showConfidence, setShowConfidence] = useState(true);

  return (
    <section className="page results-page">
      <p className="viewing-run">
        <Folder size={15} />
        Viewing run: <strong>{selectedRun.name}</strong>
      </p>
      <PageHeader
        title="Prediction Results"
        subtitle="Ab-initio eukaryotic gene prediction from FASTA sequence"
      />

      <div className="card compact-card">
        <div className="summary-grid six">
          <Metric label="Input File" value={selectedRun.fileName} />
          <Metric label="Total Bases" value={formatNumber(selectedRun.totalBases)} />
          <Metric label="Scaffolds" value={selectedRun.scaffolds} />
          <Metric label="Genes" value={selectedRun.predictedGenes} />
          <Metric label="Exons" value={selectedRun.predictedExons} />
          <Metric label="Introns" value={selectedRun.predictedIntrons} />
        </div>
      </div>

      <div className="card genome-viewer">
        <div className="viewer-controls">
          <label>
            Scaffold
            <input
              onChange={(event) => setScaffold(event.target.value)}
              value={scaffold}
            />
          </label>
          <label>
            Start
            <input
              onChange={(event) => setStart(Number(event.target.value))}
              type="number"
              value={start}
            />
          </label>
          <label>
            End
            <input
              onChange={(event) => setEnd(Number(event.target.value))}
              type="number"
              value={end}
            />
          </label>
          <button className="small-primary">Go</button>
          <div className="viewer-tools">
            <button aria-label="Zoom out">
              <Search size={15} />
            </button>
            <button aria-label="Zoom in">
              <Search size={15} />
            </button>
          </div>
        </div>

        <GenomeTrack
          predictions={predictions}
          start={start}
          end={end}
        />

        <ConfidencePanel
          end={end}
          predictions={predictions}
          showConfidence={showConfidence}
          start={start}
          toggleConfidence={() => setShowConfidence((visible) => !visible)}
        />
      </div>

      <div className="results-grid">
        <div className="card table-card">
          <div className="table-header">
            <h2>Predicted Genes</h2>
            <div className="table-search">
              <Search size={15} />
              <input placeholder="Search prediction ID" />
            </div>
          </div>
          <table>
            <thead>
              <tr>
                <th>ID</th>
                <th>Scaffold</th>
                <th>Start</th>
                <th>End</th>
                <th>Len</th>
                <th>Exons</th>
                <th>Introns</th>
              </tr>
            </thead>
            <tbody>
              {predictions.map((prediction) => (
                <tr
                  className={
                    prediction.id === selectedPredictionId ? "selected" : ""
                  }
                  key={prediction.id}
                  onClick={() => setSelectedPredictionId(prediction.id)}
                >
                  <td>{prediction.id}</td>
                  <td>{prediction.scaffold}</td>
                  <td>{prediction.start}</td>
                  <td>{prediction.end}</td>
                  <td>{prediction.end - prediction.start + 1}</td>
                  <td>{prediction.exons.length}</td>
                  <td>{prediction.introns.length}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>

        <SelectedPrediction prediction={selectedPrediction} />
      </div>

      <div className="card export-card">
        <SectionTitle icon={<Download size={20} />} title="Export Predictions" />
        <div className="export-options">
          <label>
            <input defaultChecked type="checkbox" />
            GFF3 predictions
          </label>
          <label>
            <input type="checkbox" />
            CSV coordinate table
          </label>
          <label>
            <input type="checkbox" />
            BED coordinate file
          </label>
          <label>
            <input type="checkbox" />
            Predicted nucleotide FASTA
          </label>
        </div>
        <div className="export-actions">
          <button>
            <FileCode2 size={16} />
            Export Selected Prediction
          </button>
          <button className="primary-action inline">
            <FileArchive size={16} />
            Export All Predictions
          </button>
        </div>
      </div>
    </section>
  );
}

function ConfidencePanel({
  end,
  predictions,
  showConfidence,
  start,
  toggleConfidence,
}) {
  const rangeLength = Math.max(end - start + 1, 0);
  const canShowConfidence = rangeLength > 0 && rangeLength <= 100;
  const confidenceData = useMemo(
    () =>
      Array.from({ length: rangeLength }, (_, index) => {
        const position = start + index;
        const state = predictedStateAt(position, predictions);
        return {
          confidence: mockConfidenceFor(position, state),
          position,
          state,
        };
      }),
    [rangeLength, start, predictions]
  );

  return (
    <div className="confidence-panel">
      <div className="confidence-header">
        <div>
          <h3>Per-Base Prediction Confidence</h3>
          <p>Planned Forward-Backward posterior confidence display</p>
        </div>
        <label className="toggle-row">
          <span>Show per-base confidence</span>
          <button
            aria-pressed={showConfidence}
            className={showConfidence ? "toggle on" : "toggle"}
            onClick={toggleConfidence}
            type="button"
          >
            <span />
          </button>
        </label>
      </div>

      {!showConfidence ? (
        <div className="confidence-empty">
          Enable the toggle to preview the confidence histogram.
        </div>
      ) : !canShowConfidence ? (
        <div className="confidence-empty">
          Select a range of 100 bp or less to view per-base confidence.
        </div>
      ) : (
        <>
          <div className="confidence-chart" aria-label="Per-base confidence histogram">
            <div className="confidence-y-axis">
              <span>1.0</span>
              <span>0.5</span>
              <span>0.0</span>
            </div>
            <div className="confidence-bars">
              {confidenceData.map((base) => (
                <div
                  className={`confidence-bar ${stateClass(base.state)}`}
                  key={base.position}
                  style={{ height: `${Math.max(base.confidence * 100, 4)}%` }}
                  title={`Position: ${base.position}\nState: ${base.state}\nConfidence: ${base.confidence.toFixed(2)}`}
                />
              ))}
            </div>
          </div>
          <div className="confidence-axis">
            <span>{start}</span>
            <span>Genomic position (bp)</span>
            <span>{end}</span>
          </div>
          <p className="confidence-note">
            Placeholder values are shown until the model outputs posterior
            probabilities. Confidence should represent the model's posterior
            confidence in the predicted state at each base, not biological
            accuracy or validation against a reference annotation.
          </p>
        </>
      )}
    </div>
  );
}

function GenomeTrack({ predictions, start, end }) {
  const visiblePredictions = predictions.filter(
    (prediction) => prediction.end >= start && prediction.start <= end
  );
  const range = Math.max(end - start, 1);

  function xFor(position) {
    return 5 + ((position - start) / range) * 90;
  }

  return (
    <div className="track-wrap">
      <div className="axis-labels">
        <span>{start}</span>
        <span>{Math.round(start + range * 0.33)}</span>
        <span>{Math.round(start + range * 0.66)}</span>
        <span>{end}</span>
      </div>
      <svg className="gene-track" viewBox="0 0 100 22" preserveAspectRatio="none">
        <line className="intergenic-line" x1="2" x2="98" y1="11" y2="11" />
        {visiblePredictions.map((prediction) => {
          const geneStart = Math.max(xFor(prediction.start), 2);
          const geneEnd = Math.min(xFor(prediction.end), 98);
          return (
            <g key={prediction.id}>
              <line
                className="intron-line"
                x1={geneStart}
                x2={geneEnd}
                y1="11"
                y2="11"
              />
              {prediction.exons.map(([exonStart, exonEnd]) => (
                <rect
                  className="exon-block"
                  height="5"
                  key={`${prediction.id}-${exonStart}`}
                  rx="0.35"
                  width={Math.max(xFor(exonEnd) - xFor(exonStart), 1)}
                  x={Math.max(xFor(exonStart), 2)}
                  y="8.5"
                />
              ))}
              <line
                className="start-marker"
                x1={geneStart}
                x2={geneStart}
                y1="5"
                y2="17"
              />
              <line
                className="stop-marker"
                x1={geneEnd}
                x2={geneEnd}
                y1="5"
                y2="17"
              />
            </g>
          );
        })}
      </svg>
      <div className="track-legend">
        <span>
          <i className="legend-gray" /> Intergenic
        </span>
        <span>
          <i className="legend-blue" /> Exon
        </span>
        <span>
          <i className="legend-purple" /> Intron
        </span>
        <span>
          <i className="legend-green" /> Start
        </span>
        <span>
          <i className="legend-red" /> Stop
        </span>
      </div>
    </div>
  );
}

function SelectedPrediction({ prediction }) {
  if (!prediction) {
    return (
      <div className="card selected-card">
        <h2>Selected Prediction</h2>
        <p>No prediction selected.</p>
      </div>
    );
  }

  return (
    <div className="card selected-card">
      <div className="selected-title">
        <h2>Selected Prediction</h2>
        <button aria-label="Copy sequence">
          <Copy size={15} />
        </button>
      </div>
      <dl>
        <dt>ID</dt>
        <dd>{prediction.id}</dd>
        <dt>Locus</dt>
        <dd>
          {prediction.scaffold}:{prediction.start}-{prediction.end}
        </dd>
        <dt>Length</dt>
        <dd>{prediction.end - prediction.start + 1} bp</dd>
        <dt>Exon intervals</dt>
        <dd>{formatIntervals(prediction.exons)}</dd>
        <dt>Intron intervals</dt>
        <dd>{prediction.introns.length ? formatIntervals(prediction.introns) : "None"}</dd>
        <dt>Start codon position</dt>
        <dd>{prediction.start}</dd>
        <dt>Stop codon position</dt>
        <dd>{prediction.end - 2}</dd>
      </dl>
      <h3>Sequence Preview</h3>
      <pre>{highlightSequence(prediction.sequence)}</pre>
    </div>
  );
}

function PageHeader({ title, subtitle }) {
  return (
    <header className="page-header">
      <h1>{title}</h1>
      <p>{subtitle}</p>
    </header>
  );
}

function SectionTitle({ icon, title }) {
  return (
    <div className="section-title">
      {icon}
      <h2>{title}</h2>
    </div>
  );
}

function Metric({ label, value }) {
  return (
    <div className="metric">
      <span>{label}</span>
      <strong>{value}</strong>
    </div>
  );
}

function formatBytes(bytes) {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}

function formatNumber(value) {
  if (typeof value !== "number") return value;
  return new Intl.NumberFormat("en-US").format(value);
}

function formatIntervals(intervals) {
  return intervals.map(([start, end]) => `${start}-${end}`).join(", ");
}

function highlightSequence(sequence) {
  return `${sequence.slice(0, 3)} ${sequence.slice(3, -3)} ${sequence.slice(-3)}`;
}

function predictedStateAt(position, predictions) {
  const prediction = predictions.find(
    (candidate) => position >= candidate.start && position <= candidate.end
  );
  if (!prediction) return "Intergenic";
  if (position >= prediction.start && position <= prediction.start + 2) return "Start";
  if (position >= prediction.end - 2 && position <= prediction.end) return "Stop";
  if (
    prediction.introns.some(
      ([intronStart, intronEnd]) => position >= intronStart && position <= intronEnd
    )
  ) {
    return "Intron";
  }
  return "Exon";
}

function mockConfidenceFor(position, state) {
  const wave = Math.sin(position * 0.17) * 0.08;
  const baseByState = {
    Intergenic: 0.44,
    Exon: 0.82,
    Intron: 0.72,
    Start: 0.88,
    Stop: 0.86,
  };
  return Math.min(0.98, Math.max(0.18, baseByState[state] + wave));
}

function stateClass(state) {
  return state.toLowerCase();
}

export default App;
