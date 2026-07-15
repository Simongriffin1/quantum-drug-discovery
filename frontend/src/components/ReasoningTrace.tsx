/** Agent reasoning-trace panel. */

"use client";

import type { Trace } from "@/lib/types";
import { ProvenanceBadge } from "./ProvenanceBadge";

export function ReasoningTrace({ trace }: { trace: Trace | null }) {
  return (
    <section className="panel" data-testid="reasoning-trace">
      <header className="panel-header">
        <h2>Agent trace</h2>
        <span className="panel-meta">{trace ? trace.session_id.slice(0, 8) : "—"}</span>
      </header>
      {!trace ? (
        <p className="px-3 py-8 text-sm text-muted">No agent trace.</p>
      ) : (
        <div className="max-h-[320px] space-y-1 overflow-y-auto px-3 pb-3">
          {trace.events.map((e, i) => (
            <div
              key={`${e.kind}-${i}`}
              className="border-b border-panel-border/60 py-1.5 font-mono text-[11px]"
              data-testid="trace-event"
              data-kind={e.kind}
            >
              <span className="mr-2 uppercase tracking-wide text-accent">{e.kind}</span>
              {e.tool_name ? <span className="mr-2 text-muted">{e.tool_name}</span> : null}
              <span className="text-ink/90">{e.content}</span>
            </div>
          ))}
          {trace.summary ? (
            <pre className="mt-2 whitespace-pre-wrap rounded border border-panel-border bg-panel-inset p-2 text-[11px]">
              {trace.summary}
            </pre>
          ) : null}
          <div className="pt-2">
            <ProvenanceBadge provenance={trace.provenance} />
          </div>
        </div>
      )}
    </section>
  );
}
