import {
  BarChart3,
  ChevronDown,
  ChevronLeft,
  ChevronRight,
  Copy,
  Download,
  FileArchive,
  FileCode2,
  FileText,
  FlaskConical,
  Folder,
  FolderOpen,
  Play,
  Plus,
  Search,
  UploadCloud,
} from "lucide-react";
import { useEffect, useMemo, useRef, useState } from "react";

const API_BASE = "http://localhost:5174";

function intronCount(prediction) {
  return prediction?.intronCount ?? prediction?.intron_count ?? prediction?.introns?.length ?? 0;
}
const GENES_PER_PAGE = 10;
const MAX_CONFIDENCE_BP = 500;
const MAX_REGIONAL_CONFIDENCE_BP = 5000;
const DEFAULT_VIEW_BP = 100000;
const MIN_VIEW_BP = 5;
const MAX_VIEW_BP = 500000;

async function readJsonResponse(response) {
  const text = await response.text();
  if (!text) {
    throw new Error(
      response.status === 404
        ? "API route not found. Restart with: cd frontend && npm run dev:all"
        : `Empty response from API (HTTP ${response.status}).`
    );
  }
  try {
    return JSON.parse(text);
  } catch {
    throw new Error(`Invalid API response (HTTP ${response.status}).`);
  }
}

function App() {
  const [activePage, setActivePage] = useState("run");
  const [runsExpanded, setRunsExpanded] = useState(true);
  const [runs, setRuns] = useState([]);
  const [selectedRunId, setSelectedRunId] = useState("");
  const [selectedPredictionId, setSelectedPredictionId] = useState("");
  const [uploadedFiles, setUploadedFiles] = useState([]);
  const [runStatus, setRunStatus] = useState("Not Started");
  const [elapsedTime, setElapsedTime] = useState("0s");
  const [errorMessage, setErrorMessage] = useState("");
  const fileInputRef = useRef(null);

  const selectedRun = useMemo(
    () => runs.find((run) => run.id === selectedRunId) ?? null,
    [runs, selectedRunId]
  );

  const predictions = selectedRun?.predictions ?? [];
  const selectedPrediction =
    predictions.find((prediction) => prediction.id === selectedPredictionId) ??
    predictions[0] ??
    null;

  function handleNewAnalysis() {
    setActivePage("run");
    setSelectedRunId("");
    setSelectedPredictionId("");
    setUploadedFiles([]);
    setRunStatus("Not Started");
    setElapsedTime("0s");
    setErrorMessage("");
  }

  function openSavedRun(runId) {
    const run = runs.find((candidate) => candidate.id === runId);
    setSelectedRunId(runId);
    setSelectedPredictionId(run?.predictions?.[0]?.id ?? "");
    setActivePage("results");
  }

  async function handleRunPrediction() {
    const filesToRun = uploadedFiles.filter((item) => item.checked);
    if (filesToRun.length === 0) {
      setErrorMessage("Select at least one FASTA/FNA file before running prediction.");
      return;
    }

    const started = Date.now();
    setRunStatus("Running");
    setElapsedTime("0s");
    setErrorMessage("");

    const timer = window.setInterval(() => {
      setElapsedTime(`${Math.floor((Date.now() - started) / 1000)}s`);
    }, 500);

    try {
      const completedRuns = [];

      for (const item of filesToRun) {
        const body = new FormData();
        body.append("file", item.file);

        const response = await fetch(`${API_BASE}/api/predict`, {
          body,
          method: "POST",
        });
        const payload = await readJsonResponse(response);
        if (!response.ok) {
          throw new Error(payload.detail || payload.error || `Prediction failed for ${item.name}.`);
        }
        completedRuns.push(payload);
      }

      setRuns((existing) => [...completedRuns, ...existing]);
      setSelectedRunId(completedRuns[0].id);
      setSelectedPredictionId(completedRuns[0].predictions?.[0]?.id ?? "");
      setRunStatus("Complete");
      setElapsedTime(`${Math.max(1, Math.round((Date.now() - started) / 1000))}s`);
      setActivePage("results");
    } catch (error) {
      setRunStatus("Failed");
      setErrorMessage(
        error.message === "Failed to fetch"
          ? "Cannot reach the local API (port 5174). From frontend/, run: npm run dev:all"
          : error.message
      );
    } finally {
      window.clearInterval(timer);
    }
  }

  async function handleFileChange(event) {
    const files = Array.from(event.target.files ?? []);
    if (files.length === 0) return;

    const parsedFiles = await Promise.all(
      files.map(async (file) => {
        const text = await file.text();
        const sequenceCount = (text.match(/^>/gm) ?? []).length || 1;
        const totalBases = text
          .split("\n")
          .filter((line) => !line.startsWith(">"))
          .join("")
          .replace(/[^ACGTNacgtn]/g, "").length;

        return {
          checked: true,
          file,
          id: `${file.name}-${file.size}-${file.lastModified}-${Math.random()}`,
          name: file.name,
          scaffolds: sequenceCount,
          size: formatBytes(file.size),
          totalBases,
        };
      })
    );

    setUploadedFiles((existing) => [...existing, ...parsedFiles]);
    setRunStatus("Ready");
    setElapsedTime("0s");
    setErrorMessage("");
  }

  function toggleUploadedFile(fileId) {
    setUploadedFiles((existing) =>
      existing.map((item) =>
        item.id === fileId ? { ...item, checked: !item.checked } : item
      )
    );
  }

  return (
    <div className="app-shell">
      <aside className="sidebar">
        <div className="brand">
          <div className="brand-mark">H</div>
          <div>
            <h1>HMM Gene Predictor</h1>
            <p>local backend</p>
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
              {runs.length === 0 ? (
                <p className="empty-runs">No local runs yet.</p>
              ) : (
                runs.map((run) => (
                  <button
                    className={
                      selectedRunId === run.id ? "saved-run active" : "saved-run"
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
                ))
              )}
            </div>
          )}
        </div>

        <div className="sidebar-footer">
          <a
            href="https://github.com/sansen405/gene-annotation-hidden-markov-model"
            rel="noopener noreferrer"
            target="_blank"
          >
            <FileText size={17} />
            Documentation
          </a>
        </div>
      </aside>

      <main className="content">
        {activePage === "run" ? (
          <RunPredictionPage
            elapsedTime={elapsedTime}
            errorMessage={errorMessage}
            fileInputRef={fileInputRef}
            onFileChange={handleFileChange}
            onRunPrediction={handleRunPrediction}
            onToggleUploadedFile={toggleUploadedFile}
            runStatus={runStatus}
            uploadedFiles={uploadedFiles}
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
  errorMessage,
  fileInputRef,
  onFileChange,
  onRunPrediction,
  onToggleUploadedFile,
  runStatus,
  uploadedFiles,
}) {
  const checkedFiles = uploadedFiles.filter((item) => item.checked);
  const checkedCount = checkedFiles.length;
  const totalBases = checkedFiles.reduce((sum, item) => sum + item.totalBases, 0);
  const totalScaffolds = checkedFiles.reduce((sum, item) => sum + item.scaffolds, 0);

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
            multiple
            onChange={onFileChange}
            ref={fileInputRef}
            type="file"
          />
        </div>

        <div className="file-summary multi-file-summary">
          {uploadedFiles.length === 0 ? (
            <div className="file-row">
              <div>
                <FileText size={18} />
                <strong>No file selected</strong>
              </div>
              <span>--</span>
            </div>
          ) : (
            uploadedFiles.map((item) => (
              <label className="file-row selectable-file-row" key={item.id}>
                <div>
                  <input
                    checked={item.checked}
                    onChange={() => onToggleUploadedFile(item.id)}
                    type="checkbox"
                  />
                  <FileText size={18} />
                  <strong>{item.name}</strong>
                </div>
                <span>{item.size}</span>
              </label>
            ))
          )}
          <div className="summary-grid two">
            <Metric
              label="Selected Files"
              value={uploadedFiles.length === 0 ? "--" : `${checkedCount}/${uploadedFiles.length}`}
            />
            <Metric
              label="Selected Bases"
              value={uploadedFiles.length === 0 ? "--" : formatNumber(totalBases)}
            />
          </div>
          {uploadedFiles.length > 0 && (
            <p className="file-summary-note">
              Total scaffolds across selected uploads: {formatNumber(totalScaffolds)}
            </p>
          )}
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
              <strong>Input:</strong>{" "}
              {uploadedFiles.length === 0
                ? "No file selected"
                : `${checkedCount} of ${uploadedFiles.length} files selected`}
            </p>
            {runStatus === "Running" && (
              <p className="running-copy">Training HMM and decoding sequence locally...</p>
            )}
          </div>
          <div className="elapsed">
            <span>Elapsed Time</span>
            <strong>{elapsedTime}</strong>
          </div>
        </div>

        {errorMessage && <p className="error-message">{errorMessage}</p>}

        <button
          className="primary-action"
          disabled={runStatus === "Running"}
          onClick={onRunPrediction}
        >
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
  const firstScaffold = selectedRun?.scaffolds?.[0]?.name ?? "";
  const [scaffold, setScaffold] = useState(firstScaffold);
  const [start, setStart] = useState(1);
  const [end, setEnd] = useState(DEFAULT_VIEW_BP);
  const [frameStartInput, setFrameStartInput] = useState("1");
  const [frameEndInput, setFrameEndInput] = useState(String(DEFAULT_VIEW_BP));
  const [showConfidence, setShowConfidence] = useState(true);
  const [genePage, setGenePage] = useState(1);
  const [geneSearch, setGeneSearch] = useState("");
  const [geneRangeStart, setGeneRangeStart] = useState("");
  const [geneRangeEnd, setGeneRangeEnd] = useState("");
  const [rangeMessage, setRangeMessage] = useState("");

  const activeScaffold = scaffold || firstScaffold;
  function scaffoldLengthFor(scaffoldName) {
    const match = selectedRun?.scaffolds?.find((item) => item.name === scaffoldName);
    return match?.length ?? DEFAULT_VIEW_BP;
  }

  const scaffoldLength = scaffoldLengthFor(activeScaffold);

  useEffect(() => {
    const defaultEnd = Math.min(DEFAULT_VIEW_BP, scaffoldLengthFor(firstScaffold));
    setScaffold(firstScaffold);
    setStart(1);
    setEnd(defaultEnd);
    setFrameStartInput("1");
    setFrameEndInput(String(defaultEnd));
    setGenePage(1);
    setGeneSearch("");
    setGeneRangeStart("");
    setGeneRangeEnd("");
    setRangeMessage("");
  }, [firstScaffold, selectedRun?.id]);

  function clampRangeForLength(nextStart, nextEnd, length) {
    const boundedLength = Math.max(1, length);
    const requestedStart = Number(nextStart) || 1;
    const requestedEnd = Number(nextEnd) || requestedStart;
    const minFrame = Math.min(MIN_VIEW_BP, boundedLength);
    const maxFrame = Math.min(MAX_VIEW_BP, boundedLength);
    const messages = [];

    let normalizedStart = Math.max(1, Math.min(requestedStart, boundedLength));
    let normalizedEnd = Math.max(1, Math.min(requestedEnd, boundedLength));

    if (requestedStart < 1 || requestedEnd < 1 || requestedStart > boundedLength || requestedEnd > boundedLength) {
      messages.push(`Coordinates must stay within this scaffold (1-${formatNumber(boundedLength)} bp).`);
    }

    if (normalizedEnd < normalizedStart) {
      normalizedEnd = normalizedStart;
      messages.push("Start must be before end; adjusted the frame.");
    }

    let frameLength = normalizedEnd - normalizedStart + 1;
    if (frameLength < minFrame) {
      normalizedEnd = Math.min(boundedLength, normalizedStart + minFrame - 1);
      normalizedStart = Math.max(1, normalizedEnd - minFrame + 1);
      frameLength = normalizedEnd - normalizedStart + 1;
      messages.push(`Visualization frame must be at least ${MIN_VIEW_BP} bp.`);
    }

    if (frameLength > maxFrame) {
      normalizedEnd = Math.min(boundedLength, normalizedStart + maxFrame - 1);
      normalizedStart = Math.max(1, normalizedEnd - maxFrame + 1);
      messages.push(`Visualization frame cannot exceed ${formatNumber(MAX_VIEW_BP)} bp.`);
    }

    return { message: messages[0] ?? "", start: normalizedStart, end: normalizedEnd };
  }

  function clampRange(nextStart, nextEnd) {
    return clampRangeForLength(nextStart, nextEnd, scaffoldLength);
  }

  function validateFrameRange(nextStart, nextEnd, length) {
    const trimmedStart = String(nextStart).trim();
    const trimmedEnd = String(nextEnd).trim();
    const numericStart = Number(trimmedStart);
    const numericEnd = Number(trimmedEnd);
    const boundedLength = Math.max(1, length);

    if (!trimmedStart || !trimmedEnd || !Number.isFinite(numericStart) || !Number.isFinite(numericEnd)) {
      return { message: "Enter numeric start and end positions.", valid: false };
    }
    if (numericStart < 1 || numericEnd > boundedLength || numericEnd < 1 || numericStart > boundedLength) {
      return {
        message: `Coordinates must stay within this scaffold (1-${formatNumber(boundedLength)} bp).`,
        valid: false,
      };
    }
    if (numericEnd < numericStart) {
      return { message: "Start must be before end.", valid: false };
    }

    const frameLength = numericEnd - numericStart + 1;
    if (frameLength < MIN_VIEW_BP) {
      return { message: `Visualization frame must be at least ${MIN_VIEW_BP} bp.`, valid: false };
    }
    if (frameLength > MAX_VIEW_BP) {
      return { message: `Visualization frame cannot exceed ${formatNumber(MAX_VIEW_BP)} bp.`, valid: false };
    }

    return { end: numericEnd, message: "", start: numericStart, valid: true };
  }

  function setCommittedFrame(nextStart, nextEnd, message = "") {
    setStart(nextStart);
    setEnd(nextEnd);
    setFrameStartInput(String(nextStart));
    setFrameEndInput(String(nextEnd));
    setRangeMessage(message);
  }

  function handleFrameInputChange(field, value) {
    const nextStartInput = field === "start" ? value : frameStartInput;
    const nextEndInput = field === "end" ? value : frameEndInput;
    if (field === "start") {
      setFrameStartInput(value);
    } else {
      setFrameEndInput(value);
    }

    const validation = validateFrameRange(nextStartInput, nextEndInput, scaffoldLength);
    setRangeMessage(validation.message);
    if (validation.valid) {
      setStart(validation.start);
      setEnd(validation.end);
    }
  }

  function handleTrackRangeSelect(nextStart, nextEnd) {
    const next = clampRange(nextStart, nextEnd);
    const validation = validateFrameRange(next.start, next.end, scaffoldLength);
    if (!validation.valid) {
      setRangeMessage(validation.message);
      return;
    }
    setCommittedFrame(next.start, next.end, next.message);
  }

  function handlePredictionSelect(prediction) {
    const nextScaffoldLength = scaffoldLengthFor(prediction.scaffold);
    const geneLength = Math.max(prediction.end - prediction.start + 1, 1);
    const margin = Math.ceil(geneLength * 0.25);
    const next = clampRangeForLength(
      prediction.start - margin,
      prediction.end + margin,
      nextScaffoldLength
    );

    setSelectedPredictionId(prediction.id);
    setScaffold(prediction.scaffold);
    setCommittedFrame(next.start, next.end, next.message);
  }

  const filteredPredictions = useMemo(() => {
    const query = geneSearch.trim().toLowerCase();
    const rangeStart = Number(geneRangeStart);
    const rangeEnd = Number(geneRangeEnd);
    const hasRangeStart = geneRangeStart !== "" && Number.isFinite(rangeStart);
    const hasRangeEnd = geneRangeEnd !== "" && Number.isFinite(rangeEnd);
    const minBp = hasRangeStart && hasRangeEnd ? Math.min(rangeStart, rangeEnd) : rangeStart;
    const maxBp = hasRangeStart && hasRangeEnd ? Math.max(rangeStart, rangeEnd) : rangeEnd;

    return predictions.filter(
      (prediction) => {
        const matchesQuery =
          !query ||
        prediction.id.toLowerCase().includes(query) ||
          prediction.scaffold.toLowerCase().includes(query);
        const matchesRange =
          (!hasRangeStart || prediction.end >= minBp) &&
          (!hasRangeEnd || prediction.start <= maxBp);
        return matchesQuery && matchesRange;
      }
    );
  }, [geneRangeEnd, geneRangeStart, geneSearch, predictions]);

  const totalGenePages = Math.max(1, Math.ceil(filteredPredictions.length / GENES_PER_PAGE));

  useEffect(() => {
    setGenePage((page) => Math.min(page, totalGenePages));
  }, [totalGenePages]);

  const pagedPredictions = filteredPredictions.slice(
    (genePage - 1) * GENES_PER_PAGE,
    genePage * GENES_PER_PAGE
  );
  const confidenceData = useMemo(() => {
    const allConfidence = selectedRun?.confidenceByScaffold?.[activeScaffold] ?? [];
    return allConfidence.filter(
      (base) => base.position >= start && base.position <= end
    );
  }, [activeScaffold, end, selectedRun, start]);

  if (!selectedRun) {
    return (
      <section className="page results-page">
        <PageHeader
          title="Prediction Results"
          subtitle="Run a local prediction to view results."
        />
        <div className="card empty-results">
          <Folder size={30} />
          <h2>No prediction run selected</h2>
          <p>Upload a FASTA/FNA file and run the HMM predictor to populate this page.</p>
        </div>
      </section>
    );
  }

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
          <Metric label="Total Bases" value={formatNumber(selectedRun.summary.totalBases)} />
          <Metric label="Scaffolds" value={selectedRun.summary.scaffolds} />
          <Metric label="Genes" value={selectedRun.summary.genes} />
          <Metric label="UTR Segments" value={selectedRun.summary.exons} />
          <Metric label="Total Introns" value={selectedRun.summary.introns} />
        </div>
      </div>

      <div className="card genome-viewer">
        <div className="viewer-controls">
          <label>
            Scaffold
            <select
              onChange={(event) => {
                const nextScaffold = event.target.value;
                const nextLength =
                  selectedRun.scaffolds.find((item) => item.name === nextScaffold)?.length ??
                  DEFAULT_VIEW_BP;
                const nextEnd = Math.min(DEFAULT_VIEW_BP, nextLength);
                setScaffold(nextScaffold);
                setCommittedFrame(1, nextEnd);
              }}
              value={activeScaffold}
            >
              {selectedRun.scaffolds.map((item) => (
                <option key={item.name} value={item.name}>
                  {item.name}
                </option>
              ))}
            </select>
          </label>
          <label>
            Start
            <input
              onBlur={() => setFrameStartInput(String(start))}
              onChange={(event) => handleFrameInputChange("start", event.target.value)}
              type="number"
              value={frameStartInput}
            />
          </label>
          <label>
            End
            <input
              onBlur={() => setFrameEndInput(String(end))}
              onChange={(event) => handleFrameInputChange("end", event.target.value)}
              type="number"
              value={frameEndInput}
            />
          </label>
        </div>
        {rangeMessage ? (
          <div className="track-wrap">
            <div className="frame-error">{rangeMessage}</div>
          </div>
        ) : (
          <GenomeTrack
            onRangeSelect={handleTrackRangeSelect}
            predictions={predictions.filter((item) => item.scaffold === activeScaffold)}
            start={start}
            end={end}
          />
        )}

        <ConfidencePanel
          confidenceData={confidenceData}
          end={end}
          rangeMessage={rangeMessage}
          showConfidence={showConfidence}
          start={start}
          toggleConfidence={() => setShowConfidence((visible) => !visible)}
        />
      </div>

      <div className="results-grid">
        <div className="card table-card">
          <div className="table-header">
            <h2>Predicted Genes</h2>
            <div className="table-filters">
              <label>
                <span>Search:</span>
                <span className="table-filter-input">
                  <Search size={15} />
                  <input
                    onChange={(event) => {
                      setGeneSearch(event.target.value);
                      setGenePage(1);
                    }}
                    value={geneSearch}
                  />
                </span>
              </label>
              <label className="bp-range-filter">
                <span>BP Range:</span>
                <input
                  min={1}
                  onChange={(event) => {
                    setGeneRangeStart(event.target.value);
                    setGenePage(1);
                  }}
                  type="number"
                  value={geneRangeStart}
                />
                <span>-</span>
                <input
                  min={1}
                  onChange={(event) => {
                    setGeneRangeEnd(event.target.value);
                    setGenePage(1);
                  }}
                  type="number"
                  value={geneRangeEnd}
                />
              </label>
            </div>
          </div>
          <div className="table-scroll">
            <table>
              <thead>
                <tr>
                  <th>ID</th>
                  <th>Scaffold</th>
                  <th>Strand</th>
                  <th>Start</th>
                  <th>End</th>
                  <th>Len</th>
                  <th>UTR</th>
                  <th>Introns</th>
                </tr>
              </thead>
              <tbody>
                {pagedPredictions.map((prediction) => (
                  <tr
                    className={
                      prediction.id === selectedPredictionId ? "selected" : ""
                    }
                    key={prediction.id}
                    onClick={() => handlePredictionSelect(prediction)}
                  >
                    <td>{prediction.id}</td>
                    <td>{prediction.scaffold}</td>
                    <td>{prediction.strand ?? "+"}</td>
                    <td>{prediction.start}</td>
                    <td>{prediction.end}</td>
                    <td>{prediction.end - prediction.start + 1}</td>
                    <td>{prediction.exons.length}</td>
                    <td>{intronCount(prediction)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <div className="table-pagination">
            <span>
              {filteredPredictions.length === 0
                ? "No genes"
                : `${(genePage - 1) * GENES_PER_PAGE + 1}-${Math.min(
                    genePage * GENES_PER_PAGE,
                    filteredPredictions.length
                  )} of ${filteredPredictions.length}`}
            </span>
            <div className="pagination-controls">
              <button
                aria-label="Previous page"
                disabled={genePage <= 1}
                onClick={() => setGenePage((page) => Math.max(1, page - 1))}
                type="button"
              >
                <ChevronLeft size={16} />
              </button>
              <span>
                Page {genePage} / {totalGenePages}
              </span>
              <button
                aria-label="Next page"
                disabled={genePage >= totalGenePages}
                onClick={() =>
                  setGenePage((page) => Math.min(totalGenePages, page + 1))
                }
                type="button"
              >
                <ChevronRight size={16} />
              </button>
            </div>
          </div>
        </div>

        <SelectedPrediction prediction={selectedPrediction} />
      </div>

    </section>
  );
}

function ConfidencePanel({
  confidenceData,
  end,
  rangeMessage,
  showConfidence,
  start,
  toggleConfidence,
}) {
  const rangeLength = Math.max(end - start + 1, 0);
  const canShowPerBase = rangeLength > 0 && rangeLength <= MAX_CONFIDENCE_BP;
  const canShowRegional = rangeLength > 0 && rangeLength <= MAX_REGIONAL_CONFIDENCE_BP;

  return (
    <div className="confidence-panel">
      <div className="confidence-header">
        <div>
          <h3>Prediction Confidence</h3>
          <p>
            {showConfidence
              ? "Per-base Forward-Backward posterior for the selected range"
              : "Average confidence by predicted region in the selected range"}
          </p>
        </div>
        <label className="toggle-row">
          <span>Per-base detail</span>
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

      {rangeMessage ? (
        <div className="confidence-empty error">{rangeMessage}</div>
      ) : confidenceData.length === 0 ? (
        <div className="confidence-empty">
          No confidence values are available for this range.
        </div>
      ) : showConfidence && !canShowPerBase ? (
        <div className="confidence-empty">
          {canShowRegional
            ? "Per-base detail is limited to 500 bp. Turn it off to view regional averages for this window."
            : "Per-base detail is limited to 500 bp. Regional averages are limited to 5,000 bp."}
        </div>
      ) : !showConfidence && !canShowRegional ? (
        <div className="confidence-empty">
          Regional averages are limited to 5,000 bp. Select a smaller window to view confidence.
        </div>
      ) : !showConfidence ? (
        <>
          <div className="confidence-chart grouped" aria-label="Average confidence by predicted group">
            <div className="confidence-y-axis">
              <span>1.0</span>
              <span>0.5</span>
              <span>0.0</span>
            </div>
            <div className="confidence-bars grouped-bars">
              {groupConfidenceData(confidenceData).map((group) => (
                <div className="confidence-group" key={`${group.state}-${group.start}-${group.end}`}>
                  <div className="confidence-group-bar-area">
                    <div
                      className={`confidence-bar ${stateClass(group.state)}`}
                      style={{ height: `${Math.max(group.average * 100, 4)}%` }}
                      title={`Range: ${group.start}-${group.end}\nState: ${group.state}\nAverage confidence: ${group.average.toFixed(2)}`}
                    />
                  </div>
                  <span>{group.state}</span>
                </div>
              ))}
            </div>
          </div>
          <div className="confidence-axis">
            <div className="confidence-axis-labels">
              <span className="axis-edge">{start}</span>
              <span className="axis-center">Regional averages</span>
              <span className="axis-edge">{end}</span>
            </div>
          </div>
        </>
      ) : (
        <>
          <div className="confidence-chart per-base" aria-label="Per-base confidence histogram">
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
            <div className="confidence-axis-labels">
              <span className="axis-edge">{start}</span>
              <span className="axis-center">Genomic position (bp)</span>
              <span className="axis-edge">{end}</span>
            </div>
          </div>
        </>
      )}

      <p className="confidence-note">
        Confidence is the model posterior for the predicted state at each base;
        it is not biological accuracy unless compared against a reference annotation.
      </p>
    </div>
  );
}

function GenomeTrack({ predictions, start, end, onRangeSelect }) {
  const trackRef = useRef(null);
  const [selection, setSelection] = useState(null);
  const visiblePredictions = predictions.filter(
    (prediction) => prediction.end >= start && prediction.start <= end
  );
  const span = Math.max(end - start, 0);
  const rangeLength = Math.max(end - start + 1, 1);

  function positionRatio(position) {
    if (span === 0) return 0.5;
    return (position - start) / span;
  }

  function xFor(position) {
    return 2 + positionRatio(position) * 96;
  }

  function boundaryFor(position, edge = "start") {
    const offset = edge === "end" ? 1 : 0;
    return 2 + ((position - start + offset) / rangeLength) * 96;
  }

  function intervalRect(intervalStart, intervalEnd, minWidth = 0.35) {
    if (intervalEnd < start || intervalStart > end) return null;
    const clippedStart = Math.max(intervalStart, start);
    const clippedEnd = Math.min(intervalEnd, end);
    const x = boundaryFor(clippedStart, "start");
    return {
      width: Math.max(boundaryFor(clippedEnd, "end") - x, minWidth),
      x,
    };
  }

  function bpFromClientX(clientX) {
    const rect = trackRef.current?.getBoundingClientRect();
    if (!rect?.width) return start;
    const ratio = Math.min(Math.max((clientX - rect.left) / rect.width, 0), 1);
    if (span === 0) return start;
    if (ratio <= 0) return start;
    if (ratio >= 1) return end;
    return Math.min(end, Math.max(start, Math.round(start + ratio * span)));
  }

  function handlePointerDown(event) {
    if (event.button !== 0) return;
    event.currentTarget.setPointerCapture(event.pointerId);
    const anchorBp = bpFromClientX(event.clientX);
    setSelection({ anchorBp, currentBp: anchorBp });
  }

  function handlePointerMove(event) {
    if (!selection) return;
    setSelection((current) =>
      current ? { ...current, currentBp: bpFromClientX(event.clientX) } : current
    );
  }

  function finishSelection(event) {
    if (!selection) return;
    event.currentTarget.releasePointerCapture(event.pointerId);
    const nextStart = Math.min(selection.anchorBp, selection.currentBp);
    const nextEnd = Math.max(selection.anchorBp, selection.currentBp);
    setSelection(null);
    if (nextEnd >= nextStart) {
      onRangeSelect(nextStart, nextEnd);
    }
  }

  const selectionStyle = useMemo(() => {
    if (!selection) return null;
    const leftBp = Math.min(selection.anchorBp, selection.currentBp);
    const rightBp = Math.max(selection.anchorBp, selection.currentBp);
    if (span === 0) {
      return { left: "0%", width: "100%" };
    }
    const leftPct = ((leftBp - start) / span) * 100;
    const widthPct = ((rightBp - leftBp) / span) * 100;
    return {
      left: `${Math.max(0, leftPct)}%`,
      width: `${Math.min(100 - leftPct, Math.max(0.4, widthPct))}%`,
    };
  }, [end, selection, span, start]);

  return (
    <div className="track-wrap">
      <div className="axis-labels">
        <span>{start}</span>
        <span>{span === 0 ? start : Math.round(start + span * 0.33)}</span>
        <span>{span === 0 ? start : Math.round(start + span * 0.66)}</span>
        <span>{end}</span>
      </div>
      <div className="track-interactive" ref={trackRef}>
        <svg
          className="gene-track"
          onPointerCancel={finishSelection}
          onPointerDown={handlePointerDown}
          onPointerMove={handlePointerMove}
          onPointerUp={finishSelection}
          viewBox="0 0 100 22"
          preserveAspectRatio="none"
        >
        <line className="intergenic-line" x1="2" x2="98" y1="11" y2="11" />
        {visiblePredictions.map((prediction) => {
          const showStartMarker = prediction.start >= start && prediction.start <= end;
          const showStopMarker = prediction.end >= start && prediction.end <= end;
          const geneStart = showStartMarker ? boundaryFor(prediction.start, "start") : null;
          const geneEnd = showStopMarker ? boundaryFor(prediction.end, "end") : null;
          return (
            <g key={prediction.id}>
              {prediction.introns.map((intron) => {
                const intronRect = intervalRect(intron.start, intron.end);
                if (!intronRect) return null;
                return (
                  <rect
                    className="intron-block"
                    height="2"
                    key={`${prediction.id}-intron-${intron.start}`}
                    width={intronRect.width}
                    x={intronRect.x}
                    y="10"
                  />
                );
              })}
              {prediction.exons.map((exon) => {
                const exonRect = intervalRect(exon.start, exon.end);
                if (!exonRect) return null;
                return (
                  <rect
                    className="exon-block"
                    height="5"
                    key={`${prediction.id}-${exon.start}`}
                    width={exonRect.width}
                    x={exonRect.x}
                    y="8.5"
                  />
                );
              })}
              {showStartMarker && (
                <line
                  className="start-marker"
                  x1={geneStart}
                  x2={geneStart}
                  y1="5"
                  y2="17"
                />
              )}
              {showStopMarker && (
                <line
                  className="stop-marker"
                  x1={geneEnd}
                  x2={geneEnd}
                  y1="5"
                  y2="17"
                />
              )}
            </g>
          );
        })}
        </svg>
        {selectionStyle && (
          <div className="range-selection" style={selectionStyle} />
        )}
      </div>
      <p className="track-hint">Drag on the track to select a sub-range</p>
      <div className="track-legend">
        <span>
          <i className="legend-gray" /> Intergenic
        </span>
        <span>
          <i className="legend-blue" /> UTR
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
        <button
          aria-label="Copy sequence"
          onClick={() => navigator.clipboard?.writeText(prediction.sequence || "")}
        >
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
        <dt>Strand</dt>
        <dd>{prediction.strand ?? "+"}</dd>
        <dt>Length</dt>
        <dd>{prediction.end - prediction.start + 1} bp</dd>
        <dt>UTR intervals</dt>
        <dd>{formatIntervals(prediction.exons)}</dd>
        <dt>Intron intervals</dt>
        <dd>{intronCount(prediction) || "None"}</dd>
        <dt>Start codon position</dt>
        <dd>{(prediction.strand ?? "+") === "-" ? prediction.end - 2 : prediction.start}</dd>
        <dt>Stop codon position</dt>
        <dd>{(prediction.strand ?? "+") === "-" ? prediction.start : prediction.end - 2}</dd>
      </dl>
      <h3>Sequence Preview</h3>
      <pre>{highlightSequence(prediction.sequence || "")}</pre>
    </div>
  );
}

function ExportCard({ runId }) {
  const formats = [
    ["gff3", "GFF3 predictions"],
    ["csv", "CSV coordinate table"],
    ["bed", "BED coordinate file"],
    ["fasta", "Predicted nucleotide FASTA"],
  ];

  return (
    <div className="card export-card">
      <SectionTitle icon={<Download size={20} />} title="Export Predictions" />
      <div className="export-options">
        {formats.map(([format, label]) => (
          <a
            href={`${API_BASE}/api/export/${runId}/${format}`}
            key={format}
          >
            {label}
          </a>
        ))}
      </div>
      <div className="export-actions">
        <a href={`${API_BASE}/api/export/${runId}/gff3`}>
          <FileCode2 size={16} />
          Export GFF3
        </a>
        <a className="primary-action inline" href={`${API_BASE}/api/export/${runId}/fasta`}>
          <FileArchive size={16} />
          Export FASTA
        </a>
      </div>
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
  return intervals.map((item) => `${item.start}-${item.end}`).join(", ");
}

function highlightSequence(sequence) {
  if (sequence.length <= 6) return sequence;
  return (
    <>
      <span className="start-codon">{sequence.slice(0, 3)}</span>
      {" "}
      {sequence.slice(3, -3)}
      {" "}
      <span className="stop-codon">{sequence.slice(-3)}</span>
    </>
  );
}

function stateClass(state) {
  return state.toLowerCase();
}

function groupConfidenceData(confidenceData) {
  const groups = [];

  for (const base of confidenceData) {
    const previous = groups[groups.length - 1];
    if (previous && previous.state === base.state && previous.end + 1 === base.position) {
      previous.end = base.position;
      previous.total += base.confidence;
      previous.count++;
      previous.average = previous.total / previous.count;
      continue;
    }

    groups.push({
      average: base.confidence,
      count: 1,
      end: base.position,
      start: base.position,
      state: base.state,
      total: base.confidence,
    });
  }

  return groups;
}

export default App;
