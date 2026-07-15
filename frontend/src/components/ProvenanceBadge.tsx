/** Compact provenance strip — every panel shows traceability. */

import type { Provenance } from "@/lib/types";

export function ProvenanceBadge({ provenance }: { provenance: Provenance }) {
  const sha = provenance.git_sha ? provenance.git_sha.slice(0, 8) : "uncommitted";
  const tools = Object.entries(provenance.tool_versions)
    .slice(0, 4)
    .map(([k, v]) => `${k}@${v}`)
    .join(" · ");
  return (
    <div
      className="font-mono text-[10px] leading-relaxed text-muted"
      data-testid="provenance-badge"
    >
      <span>data={provenance.data_version ?? "—"}</span>
      <span className="mx-1.5 text-panel-border">|</span>
      <span>git={sha}</span>
      {tools ? (
        <>
          <span className="mx-1.5 text-panel-border">|</span>
          <span>{tools}</span>
        </>
      ) : null}
    </div>
  );
}
