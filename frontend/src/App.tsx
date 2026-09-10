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
      {/* Nav */}
      <nav style={{ maxWidth: 960, margin: "0 auto", width: "100%", boxSizing: "border-box", padding: "24px 24px 0", display: "flex", alignItems: "center", gap: 10 }}>
        <span style={{ display: "inline-flex", alignItems: "center", justifyContent: "center", width: 34, height: 34, borderRadius: 9, background: "linear-gradient(135deg, #38bdf8, #b45309)", color: "white", fontWeight: 700, fontSize: 13 }}>CS</span>
        <span style={{ fontFamily: "'Fraunces', serif", fontWeight: 600, fontSize: 18 }}>Cross-Species Atlas</span>
      </nav>

      <section style={{ maxWidth: 960, margin: "0 auto", width: "100%", boxSizing: "border-box", padding: "32px 24px 40px", display: "grid", gridTemplateColumns: "1.1fr 0.9fr", gap: 40, alignItems: "center" }}>
        <div>
          <div style={{ textTransform: "uppercase", letterSpacing: "0.14em", fontSize: 12, fontWeight: 700, color: "#0284c7", marginBottom: 14 }}>Evolution-aware gene function</div>
          <h1 style={{ fontFamily: "'Fraunces', serif", fontWeight: 600, fontSize: 32, lineHeight: 1.15, margin: "0 0 18px" }}>
            Borrow what evolution <em style={{ fontStyle: "italic", color: "#0284c7" }}>already knows.</em>
          </h1>
          <p style={{ color: "#475569", fontSize: 16, lineHeight: 1.6, maxWidth: 480, margin: 0 }}>
            Cross-Species Atlas predicts a gene's likely function by learning from well-studied genes in other
            species — tracing the evolutionary relationships that connect them.
          </p>
        </div>
        <figure style={{ margin: 0 }}>
          <img
            src="/hero.png"
            alt="Illustration of a cross-species orthology network — interconnected gene nodes from multiple species linked by orthology edges, with DNA helix and cellular motifs representing evolution-aware functional annotation transfer via graph learning"
            style={{ width: "100%", aspectRatio: "1 / 1", objectFit: "cover", borderRadius: 18, border: "1px solid #e2e8f0", boxShadow: "0 24px 50px rgba(15,23,42,0.15)", display: "block" }}
          />
        </figure>
      </section>

      <Banner ready={ready} />

      <main id="main-content" style={{ flex: 1, maxWidth: 960, width: "100%", margin: "0 auto", padding: "2rem 1.5rem" }}>
        {/* Info cards */}
        <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(260px, 1fr))", gap: "1rem", marginBottom: "2rem" }}>
          <div style={cardStyle}>
            <h3 style={cardTitleStyle}>Built on real orthology data</h3>
            <p style={cardTextStyle}>
              Genes across species are linked by their evolutionary relationships, weighted by how confidently
              they're related.
            </p>
          </div>
          <div style={cardStyle}>
            <h3 style={cardTitleStyle}>Tested rigorously</h3>
            <p style={cardTextStyle}>
              A learned graph model is compared against a strong baseline, evaluated on species it never trained on.
            </p>
          </div>
          <div style={cardStyle}>
            <h3 style={cardTitleStyle}>Starting with cell cycle</h3>
            <p style={cardTextStyle}>
              The initial focus is <strong>cell cycle</strong> function — a well-studied biological process with
              rich reference annotation.
            </p>
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
        Built on real evolutionary orthology and gene function annotation data.
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
