import { useState } from "react";

type SearchResult = {
  gene_id: string;
  species: string;
  label: number;
};

export default function GeneSearch({
  onSelect,
  apiBase,
}: {
  onSelect: (geneId: string) => void;
  apiBase: string;
}) {
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<SearchResult[]>([]);
  const [loading, setLoading] = useState(false);
  const [message, setMessage] = useState<string | null>(null);

  async function doSearch() {
    if (!query.trim()) return;
    setLoading(true);
    setMessage(null);
    try {
      const res = await fetch(`${apiBase}/genes/search?q=${encodeURIComponent(query.trim())}`);
      const data = await res.json();
      setResults(data.results || []);
      if (data.message) setMessage(data.message);
      if ((data.results || []).length === 0 && !data.message) {
        setMessage("No matching genes found.");
      }
    } catch (e) {
      setMessage("Search failed — backend unreachable.");
    } finally {
      setLoading(false);
    }
  }

  return (
    <div
      style={{
        background: "white",
        border: "1px solid #e2e8f0",
        borderRadius: 12,
        padding: "1.25rem",
        boxShadow: "0 1px 3px rgba(0,0,0,0.06)",
      }}
    >
      <h2 style={{ fontSize: "1rem", fontWeight: 600, marginBottom: 12 }}>
        Target-Species Gene Search
      </h2>
      <div style={{ display: "flex", gap: 8 }}>
        <input
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && doSearch()}
          placeholder="e.g. TARGET_GENE_0001 or SOURCE_GENE_"
          style={{
            flex: 1,
            padding: "0.6rem 0.85rem",
            border: "1px solid #cbd5e1",
            borderRadius: 8,
            fontSize: "0.9rem",
            outline: "none",
          }}
        />
        <button
          onClick={doSearch}
          disabled={loading}
          style={{
            padding: "0.6rem 1.1rem",
            background: loading ? "#94a3b8" : "#0f172a",
            color: "white",
            border: "none",
            borderRadius: 8,
            fontWeight: 600,
            cursor: loading ? "not-allowed" : "pointer",
            fontSize: "0.9rem",
          }}
        >
          {loading ? "Searching..." : "Search"}
        </button>
      </div>

      {message && (
        <p style={{ marginTop: 10, fontSize: "0.82rem", color: "#64748b" }}>{message}</p>
      )}

      {results.length > 0 && (
        <ul style={{ listStyle: "none", marginTop: 14, display: "flex", flexDirection: "column", gap: 6 }}>
          {results.map((r) => (
            <li
              key={r.gene_id}
              onClick={() => onSelect(r.gene_id)}
              style={{
                display: "flex",
                justifyContent: "space-between",
                alignItems: "center",
                padding: "0.55rem 0.75rem",
                border: "1px solid #e2e8f0",
                borderRadius: 8,
                cursor: "pointer",
                background: "#f8fafc",
                fontSize: "0.85rem",
              }}
            >
              <span>
                <strong>{r.gene_id}</strong>{" "}
                <span
                  style={{
                    display: "inline-block",
                    fontSize: "0.7rem",
                    background: "#e0f2fe",
                    color: "#0c4a6e",
                    padding: "1px 6px",
                    borderRadius: 999,
                    marginLeft: 6,
                  }}
                >
                  taxon:{r.species}
                </span>
              </span>
              <span style={{ color: "#64748b" }} onClick={(e) => { e.stopPropagation(); onSelect(r.gene_id); }}>
                View prediction &rarr;
              </span>
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}
