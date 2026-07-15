/** Calibration reliability plot (ECE bins). */

"use client";

import type { Calibration } from "@/lib/types";
import { ProvenanceBadge } from "./ProvenanceBadge";

export function CalibrationPlot({ calibration }: { calibration: Calibration | null }) {
  return (
    <section className="panel" data-testid="calibration-plot">
      <header className="panel-header">
        <h2>Calibration</h2>
        <span className="panel-meta">
          {calibration
            ? `ECE ${calibration.ece.toFixed(3)} (gate < ${calibration.ece_threshold})`
            : "—"}
        </span>
      </header>
      {!calibration ? (
        <p className="px-3 py-8 text-sm text-muted">No calibration artifact.</p>
      ) : (
        <div className="space-y-3 px-3 pb-3">
          <div className="grid grid-cols-3 gap-2 font-mono text-[11px]">
            <Metric label="n" value={String(calibration.n)} />
            <Metric
              label="coverage"
              value={`${(calibration.empirical_coverage * 100).toFixed(1)}%`}
            />
            <Metric
              label="gate"
              value={calibration.passed ? "PASS" : "FAIL"}
              tone={calibration.passed ? "ok" : "bad"}
            />
          </div>
          <svg viewBox="0 0 320 140" className="h-[140px] w-full" data-testid="calibration-svg">
            <line x1="40" y1="120" x2="300" y2="20" stroke="#94a3b0" strokeDasharray="4 3" />
            {calibration.reliability_bins.map((b) => {
              const x = 40 + b.predicted_coverage * 260;
              const y = 120 - b.empirical_coverage * 100;
              return (
                <circle
                  key={b.bin_index}
                  cx={x}
                  cy={y}
                  r={Math.max(4, Math.min(10, b.n))}
                  fill="#1a5f4a"
                  opacity={0.85}
                />
              );
            })}
            <text x="40" y="136" className="fill-muted text-[10px]">
              predicted
            </text>
            <text x="4" y="70" className="fill-muted text-[10px]" transform="rotate(-90 4 70)">
              empirical
            </text>
          </svg>
          {calibration.notes ? (
            <p className="font-mono text-[10px] text-muted">{calibration.notes}</p>
          ) : null}
          <ProvenanceBadge provenance={calibration.provenance} />
        </div>
      )}
    </section>
  );
}

function Metric({
  label,
  value,
  tone,
}: {
  label: string;
  value: string;
  tone?: "ok" | "bad";
}) {
  const color =
    tone === "ok" ? "text-accent" : tone === "bad" ? "text-red-800" : "text-ink";
  return (
    <div className="rounded border border-panel-border bg-panel-inset px-2 py-1.5">
      <div className="text-muted">{label}</div>
      <div className={`text-sm ${color}`}>{value}</div>
    </div>
  );
}
