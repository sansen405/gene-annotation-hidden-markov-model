# ── Stage 1: Build C++ binaries ──────────────────────────────────────────────
FROM ubuntu:22.04 AS cpp-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        clang \
        nlohmann-json3-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY src/ src/
COPY train_model.cpp annotate.cpp ./

# train_model: trains from a genome profile JSON + raw data → saves model JSON
RUN clang++ -std=c++17 -O2 -Isrc \
        train_model.cpp \
        src/model/Transition_Model.cpp \
        src/model/Emission_Model.cpp \
        src/model/Model_IO.cpp \
        src/genome_profiles/Genome_Profile.cpp \
        src/parsers/FNA_Parser.cpp \
        src/parsers/GFF_Parser.cpp \
        -o /usr/local/bin/train_model

# annotate: loads a saved model JSON + FASTA → outputs gene annotation JSON
RUN clang++ -std=c++17 -O2 -Isrc \
        annotate.cpp \
        src/decoding/Viterbi.cpp \
        src/model/Emission_Model.cpp \
        src/model/Model_IO.cpp \
        src/genome_profiles/Genome_Profile.cpp \
        src/parsers/FNA_Parser.cpp \
        -o /usr/local/bin/annotate

# ── Stage 2: Final API image ──────────────────────────────────────────────────
FROM python:3.11-slim

# Runtime C++ stdlib needed by the compiled binaries
RUN apt-get update && apt-get install -y --no-install-recommends \
        libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

# Copy compiled binaries from build stage
COPY --from=cpp-builder /usr/local/bin/train_model /usr/local/bin/train_model
COPY --from=cpp-builder /usr/local/bin/annotate     /usr/local/bin/annotate

WORKDIR /app

# Install Python dependencies
COPY api/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy FastAPI app
COPY api/main.py .

# Copy pre-trained model JSON files.
# Run scripts/train_all_models.sh locally first to populate models/.
COPY models/ models/

EXPOSE 8000

ENV MODELS_DIR=/app/models
ENV ANNOTATE_BIN=/usr/local/bin/annotate

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
