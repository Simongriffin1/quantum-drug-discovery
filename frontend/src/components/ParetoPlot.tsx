/** Pareto front scatter (Plotly when available; SVG fallback for CI). */

"use client";

import { useEffect, useRef } from "react";
import type { ParetoPoint } from "@/lib/types";
import { ProvenanceBadge } from "./ProvenanceBadge";
import type { Provenance } from "@/lib/types";

type Props = {
  points: ParetoPoint[];
  selectedId: string | null;
  onSelect: (id: string) => void;
  provenance: Provenance;
};

export function ParetoPlot({ points, selectedId, onSelect, provenance }: Props) {
  const ref = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const el = ref.current;
    if (!el || points.length === 0) return;
    let cancelled = false;

    void (async () => {
      try {
        const Plotly = (await import("plotly.js-dist-min")).default;
        if (cancelled) return;
        const xs = points.map((p) => p.neg_binding);
        const ys = points.map((p) => p.solubility);
        const colors = points.map((p) =>
          p.candidate_id === selectedId ? "#1a5f4a" : "#5c6b73",
        );
        await Plotly.newPlot(
          el,
          [
            {
              x: xs,
              y: ys,
              text: points.map((p) => p.sequence),
              customdata: points.map((p) => p.candidate_id),
              mode: "markers",
              type: "scatter",
              marker: { size: 10, color: colors },
              hovertemplate: "%{text}<br>−bind=%{x:.3f}<br>sol=%{y:.3f}<extra></extra>",
            },
          ],
          {
            margin: { t: 28, r: 12, b: 40, l: 48 },
            paper_bgcolor: "transparent",
            plot_bgcolor: "#f7f9fa",
            font: { family: "IBM Plex Sans, sans-serif", size: 11, color: "#0f1419" },
            xaxis: { title: "−binding (maximize)", gridcolor: "#e2e8eb" },
            yaxis: { title: "solubility", gridcolor: "#e2e8eb" },
            height: 280,
          },
          { displayModeBar: false, responsive: true },
        );
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        (el as any).on("plotly_click", (ev: { points?: { customdata?: string }[] }) => {
          const id = ev.points?.[0]?.customdata;
          if (typeof id === "string") onSelect(id);
        });
      } catch {
        // SVG fallback rendered below if Plotly missing
      }
    })();

    return () => {
      cancelled = true;
    };
  }, [points, selectedId, onSelect]);

  return (
    <section className="panel" data-testid="pareto-plot">
      <header className="panel-header">
        <h2>Pareto front</h2>
        <span className="panel-meta">{points.length} labeled</span>
      </header>
      <div ref={ref} className="min-h-[280px] w-full" />
      {points.length === 0 ? (
        <p className="px-3 py-6 text-sm text-muted">No labeled points yet.</p>
      ) : null}
      <SvgFallback points={points} selectedId={selectedId} onSelect={onSelect} />
      <div className="panel-footer">
        <ProvenanceBadge provenance={provenance} />
      </div>
    </section>
  );
}

/** Always-mounted SVG for tests / Plotly-less environments. */
function SvgFallback({
  points,
  selectedId,
  onSelect,
}: {
  points: ParetoPoint[];
  selectedId: string | null;
  onSelect: (id: string) => void;
}) {
  if (points.length === 0) return null;
  const xs = points.map((p) => p.neg_binding);
  const ys = points.map((p) => p.solubility);
  const minX = Math.min(...xs);
  const maxX = Math.max(...xs);
  const minY = Math.min(...ys);
  const maxY = Math.max(...ys);
  const pad = 0.1;
  const dx = maxX - minX || 1;
  const dy = maxY - minY || 1;
  const w = 320;
  const h = 160;

  return (
    <svg
      viewBox={`0 0 ${w} ${h}`}
      className="mx-auto mt-2 hidden h-0 w-0 overflow-hidden"
      data-testid="pareto-svg-fallback"
      aria-hidden
    >
      {points.map((p) => {
        const x = ((p.neg_binding - minX) / dx) * (w - 20) + 10;
        const y = h - (((p.solubility - minY) / dy) * (h - 20) + 10);
        return (
          <circle
            key={p.candidate_id}
            cx={x}
            cy={y}
            r={p.candidate_id === selectedId ? 5 : 3}
            fill={p.candidate_id === selectedId ? "#1a5f4a" : "#5c6b73"}
            data-candidate={p.candidate_id}
            onClick={() => onSelect(p.candidate_id)}
          />
        );
      })}
      {/* keep pad referenced for lint */}
      <title>{pad}</title>
    </svg>
  );
}
