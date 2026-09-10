type ReadyState = {
  ready: boolean;
  approved: boolean;
  revision: string | null;
  message: string;
} | null;

export default function Banner({ ready }: { ready: ReadyState }) {
  if (ready === null) {
    return (
      <div role="status" aria-live="polite" style={{ background: "#f1f5f9", color: "#64748b", textAlign: "center", padding: "0.6rem", fontSize: "0.85rem" }}>
        Checking model status...
      </div>
    );
  }

  if (!ready.ready) {
    return (
      <div
        role="alert"
        aria-live="assertive"
        style={{
          background: "#fef3c7",
          borderBottom: "1px solid #fcd34d",
          color: "#92400e",
          textAlign: "center",
          padding: "0.75rem 1rem",
          fontSize: "0.85rem",
        }}
      >
        <strong>Predictions aren't available yet.</strong>{" "}
        <span style={{ opacity: 0.85 }}>Our team is finishing validation before enabling live results.</span>
      </div>
    );
  }

  return (
    <div
      role="status"
      aria-live="polite"
      style={{
        background: "#dcfce7",
        borderBottom: "1px solid #86efac",
        color: "#166534",
        textAlign: "center",
        padding: "0.6rem 1rem",
        fontSize: "0.85rem",
      }}
    >
      Model ready — revision <code>{ready.revision}</code> — predictions enabled.
    </div>
  );
}
