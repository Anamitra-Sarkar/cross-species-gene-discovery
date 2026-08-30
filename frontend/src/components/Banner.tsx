type ReadyState = {
  ready: boolean;
  approved: boolean;
  revision: string | null;
  message: string;
} | null;

export default function Banner({ ready }: { ready: ReadyState }) {
  if (ready === null) {
    return (
      <div style={{ background: "#f1f5f9", color: "#64748b", textAlign: "center", padding: "0.6rem", fontSize: "0.85rem" }}>
        Checking model status...
      </div>
    );
  }

  if (!ready.ready) {
    return (
      <div
        style={{
          background: "#fef3c7",
          borderBottom: "1px solid #fcd34d",
          color: "#92400e",
          textAlign: "center",
          padding: "0.75rem 1rem",
          fontSize: "0.85rem",
        }}
      >
        <strong>Model not yet released</strong> — predictions unavailable.{" "}
        <span style={{ opacity: 0.85 }}>{ready.message}</span>
        <span style={{ display: "block", fontSize: "0.78rem", marginTop: 2, opacity: 0.7 }}>
          The backend release gate is closed. Set <code>MODEL_RELEASE_APPROVED=true</code> and{" "}
          <code>APPROVED_ARTIFACT_REVISION</code> to enable predictions.
        </span>
      </div>
    );
  }

  return (
    <div
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
