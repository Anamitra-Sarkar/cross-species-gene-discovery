import { useEffect, useState } from "react";
import GeneSearch from "./components/GeneSearch";
import GeneDetail from "./components/GeneDetail";
import Banner from "./components/Banner";

type ReadyState = {
  ready: boolean;
  approved: boolean;
  revision: string | null;
  message: string;
};

export default function App() {
  const [ready, setReady] = useState<ReadyState | null>(null);
  const [selectedGene, setSelectedGene] = useState<string | null>(null);

  const apiBase = import.meta.env.VITE_API_URL || import.meta.env.VITE_API_BASE || "";

  useEffect(() => {
    fetch(`${apiBase}/ready`)
      .then((r) => r.json())
      .then(setReady)
      .catch(() =>
        setReady({ ready: false, approved: false, revision: null, message: "Backend unreachable" })
      );
  }, [apiBase]);

  return (
    <div style={{ minHeight: "100vh", display: "flex", flexDirection: "column" }}>
      {/* Header */}
      <header
        style={{
          background: "linear-gradient(135deg, #0f172a 0%, #1e3a5f 100%)",
          color: "white",
          padding: "1.5rem 2rem",
          borderBottom: "3px solid #38bdf8",
        }}
      >
        <div style={{ maxWidth: 960, margin: "0 auto" }}>
          <h1 style={{ fontSize: "1.5rem", fontWeight: 700, letterSpacing: "-0.02em" }}>
            Cross-Species Gene Function Discovery
          </h1>
          <p style={{ opacity: 0.8, fontSize: "0.9rem", marginTop: 4 }}>
            Evolution-aware graph learning — orthology-based functional annotation transfer
            via random walk &amp; GNN over cross-species orthology graphs.
          </p>
        </div>
      </header>

      <Banner ready={ready} />

      <main style={{ flex: 1, maxWidth: 960, width: "100%", margin: "0 auto", padding: "2rem 1.5rem" }}>
        {/* Info cards */}
        <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(280px, 1fr))", gap: "1rem", marginBottom: "2rem" }}>
          <div style={cardStyle}>
            <h3 style={cardTitleStyle}>Orthology Graph</h3>
            <p style={cardTextStyle}>
              Nodes are genes from 2+ species; edges connect orthologs from OrthoDB / eggNOG,
              weighted by orthology confidence.
            </p>
            <span style={badgeStyle}>OrthoDB + GO / GOA</span>
          </div>
          <div style={cardStyle}>
            <h3 style={cardTitleStyle}>Two Models</h3>
            <p style={cardTextStyle}>
              <strong>RWR</strong> (random walk with restart) baseline +{" "}
              <strong>GraphSAGE/GCN</strong> learned model. Cross-species holdout
              evaluation (source-only training).
            </p>
            <span style={badgeStyle}>RWR &amp; GNN</span>
          </div>
          <div style={cardStyle}>
            <h3 style={cardTitleStyle}>Target: GO:0007049</h3>
            <p style={cardTextStyle}>
              Default target is <strong>cell cycle</strong> (GO:0007049) — a well-studied
              Biological Process term with extensive yeast annotation.
            </p>
            <span style={badgeStyle}>GO Biological Process</span>
          </div>
        </div>

        <GeneSearch onSelect={setSelectedGene} apiBase={apiBase} />

        {selectedGene && (
          <GeneDetail geneId={selectedGene} apiBase={apiBase} ready={ready?.ready ?? false} />
        )}

        {!selectedGene && (
          <div
            style={{
              textAlign: "center",
              padding: "3rem 1rem",
              color: "#64748b",
              border: "1px dashed #cbd5e1",
              borderRadius: 12,
              background: "white",
              marginTop: "1.5rem",
            }}
          >
            <p style={{ fontSize: "1.05rem" }}>Search for a target-species gene above to view its predicted function.</p>
            <p style={{ fontSize: "0.85rem", marginTop: 8 }}>
              Predictions show GO-term scores and which source-species orthologs drove the prediction.
            </p>
          </div>
        )}
      </main>

      <footer
        style={{
          borderTop: "1px solid #e2e8f0",
          padding: "1rem 2rem",
          textAlign: "center",
          fontSize: "0.8rem",
          color: "#94a3b8",
          background: "white",
        }}
      >
        Data sources: OrthoDB (Kuznetsov et al. 2023), Gene Ontology / GOA (Gene Ontology Consortium 2023) — see{" "}
        <code>docs/data_sources.md</code> for endpoints &amp; citations. Test fixtures are explicitly synthetic.
      </footer>
    </div>
  );
}

const cardStyle: React.CSSProperties = {
  background: "white",
  border: "1px solid #e2e8f0",
  borderRadius: 12,
  padding: "1.25rem",
  boxShadow: "0 1px 3px rgba(0,0,0,0.06)",
};
const cardTitleStyle: React.CSSProperties = { fontSize: "0.95rem", fontWeight: 600, color: "#0f172a", marginBottom: 6 };
const cardTextStyle: React.CSSProperties = { fontSize: "0.85rem", color: "#475569", lineHeight: 1.5 };
const badgeStyle: React.CSSProperties = {
  display: "inline-block",
  marginTop: 10,
  fontSize: "0.7rem",
  fontWeight: 600,
  letterSpacing: "0.04em",
  textTransform: "uppercase",
  background: "#f1f5f9",
  color: "#334155",
  padding: "2px 8px",
  borderRadius: 999,
  border: "1px solid #e2e8f0",
};
