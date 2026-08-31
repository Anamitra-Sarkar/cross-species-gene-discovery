import { useEffect, useState } from "react";

type OrthologInfo = {
  gene_id: string;
  species: string;
  orthology_confidence: number;
  source_positive: boolean;
};

type Prediction = {
  gene_id: string;
  species: string;
  predicted_score: number;
  go_term: string;
  go_term_name: string;
  explanation: { method: string; source_orthologs: OrthologInfo[]; num_source_positives: number };
  model_revision: string | null;
};

export default function GeneDetail({
  geneId,
  apiBase,
  ready,
}: {
  geneId: string;
  apiBase: string;
  ready: boolean;
}) {
  const [data, setData] = useState<Prediction | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    setLoading(true);
    setError(null);
    setData(null);
    fetch(`${apiBase}/genes/${encodeURIComponent(geneId)}/predict`)
      .then(async (r) => {
        if (!r.ok) {
          const body = await r.json().catch(() => ({ detail: r.statusText }));
          throw new Error(body.detail?.message || body.detail || r.statusText);
        }
        return r.json();
      })
      .then(setData)
      .catch((e) => setError(e.message))
      .finally(() => setLoading(false));
  }, [geneId, apiBase]);

  if (loading) {
    return (
      <div style={cardStyle} aria-busy="true" aria-live="polite">
        <p style={{ color: "#64748b" }}>Loading prediction for {geneId}...</p>
      </div>
    );
  }

  if (error) {
    const isNotReleased = error.toLowerCase().includes("not yet released") || !ready;
    return (
      <div role="alert" aria-live="assertive" style={{ ...cardStyle, borderColor: isNotReleased ? "#fcd34d" : "#fecaca", background: isNotReleased ? "#fffbeb" : "#fef2f2" }}>
        <h3 style={{ color: isNotReleased ? "#92400e" : "#991b1b" }}>
          {isNotReleased ? "Model not yet released — prediction unavailable" : "Error"}
        </h3>
        <p style={{ fontSize: "0.85rem", color: "#475569", marginTop: 6 }}>{error}</p>
        <p style={{ fontSize: "0.82rem", color: "#64748b", marginTop: 6 }}>
          Gene: <code>{geneId}</code>
        </p>
      </div>
    );
  }

  if (!data) return null;

  const scoreNum = Number(data.predicted_score);
  const safeScore = Number.isFinite(scoreNum) ? Math.max(0, Math.min(1, scoreNum)) : 0;
  const pct = (safeScore * 100).toFixed(1);

  return (
    <div style={{ ...cardStyle, marginTop: "1.5rem" }} aria-live="polite">
      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "flex-start", flexWrap: "wrap", gap: 12 }}>
        <div>
          <h2 style={{ fontSize: "1.05rem", fontWeight: 700 }}>
            {data.gene_id}
            <span
              style={{
                marginLeft: 8,
                fontSize: "0.7rem",
                background: "#e0f2fe",
                color: "#0c4a6e",
                padding: "2px 7px",
                borderRadius: 999,
              }}
            >
              taxon:{data.species}
            </span>
          </h2>
          <p style={{ fontSize: "0.85rem", color: "#475569", marginTop: 4 }}>
            Predicted function: <strong>{data.go_term}</strong> — {data.go_term_name}
          </p>
        </div>
        <div style={{ textAlign: "right" }}>
          <div
            aria-label={`Predicted score ${pct} percent`}
            style={{
              display: "inline-block",
              background: safeScore > 0.5 ? "#dcfce7" : "#f1f5f9",
              color: safeScore > 0.5 ? "#166534" : "#475569",
              padding: "0.35rem 0.7rem",
              borderRadius: 8,
              fontWeight: 700,
              fontSize: "1.1rem",
            }}
          >
            {pct}%
          </div>
          <p style={{ fontSize: "0.72rem", color: "#94a3b8", marginTop: 2 }}>RWR score (0–100%)</p>
        </div>
      </div>

      {/* Explanation */}
      <div style={{ marginTop: "1rem", borderTop: "1px solid #e2e8f0", paddingTop: "1rem" }}>
        <h3 style={{ fontSize: "0.85rem", fontWeight: 600, color: "#334155", marginBottom: 8 }}>
          Which source-species orthologs drove this prediction?
        </h3>
        <p style={{ fontSize: "0.78rem", color: "#64748b", marginBottom: 10 }}>
          Method: {data.explanation.method}. {data.explanation.num_source_positives} source-positive genes used as seeds.
          Ranked by orthology confidence; highlighted rows are positively annotated for this GO term in the source species.
        </p>
        {data.explanation.source_orthologs.length === 0 ? (
          <p style={{ fontSize: "0.82rem", color: "#94a3b8" }}>No direct orthologs found for this gene in the graph.</p>
        ) : (
          <div style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            {data.explanation.source_orthologs.map((o) => (
              <div
                key={o.gene_id}
                style={{
                  display: "flex",
                  justifyContent: "space-between",
                  alignItems: "center",
                  padding: "0.5rem 0.7rem",
                  borderRadius: 8,
                  border: `1px solid ${o.source_positive ? "#86efac" : "#e2e8f0"}`,
                  background: o.source_positive ? "#f0fdf4" : "white",
                  fontSize: "0.82rem",
                }}
              >
                <span>
                  <strong>{o.gene_id}</strong>{" "}
                  <span style={{ fontSize: "0.7rem", color: "#64748b" }}>taxon:{o.species}</span>{" "}
                  {o.source_positive && (
                    <span
                      style={{
                        marginLeft: 6,
                        fontSize: "0.65rem",
                        background: "#22c55e",
                        color: "white",
                        padding: "1px 5px",
                        borderRadius: 999,
                        fontWeight: 600,
                      }}
                    >
                      SOURCE POSITIVE
                    </span>
                  )}
                </span>
                <span style={{ color: "#475569", fontWeight: 600 }}>
                  conf {o.orthology_confidence.toFixed(3)}
                </span>
              </div>
            ))}
          </div>
        )}
      </div>

      <p style={{ fontSize: "0.72rem", color: "#94a3b8", marginTop: 12 }}>
        Model revision: <code>{data.model_revision ?? "unknown"}</code> • Cross-species transfer via orthology — function may diverge
        between orthologs; prediction is probabilistic.
      </p>
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
