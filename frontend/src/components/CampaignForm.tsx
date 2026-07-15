/** Campaign start form + iteration table (interactive bench surface). */

"use client";

import type { Campaign } from "@/lib/types";
import { ProvenanceBadge } from "./ProvenanceBadge";

type Props = {
  onStart: (goal: string, seed: number, target: number) => void;
  busy: boolean;
  error: string | null;
  campaign: Campaign | null;
};

export function CampaignForm({ onStart, busy, error, campaign }: Props) {
  return (
    <section className="panel" data-testid="campaign-form">
      <header className="panel-header">
        <h2>Campaign</h2>
        <span className="panel-meta">simulation mode</span>
      </header>
      <form
        className="grid gap-3 px-3 pb-3 md:grid-cols-[1fr_auto_auto_auto]"
        onSubmit={(e) => {
          e.preventDefault();
          const fd = new FormData(e.currentTarget);
          const goal = String(fd.get("goal") || "");
          const seed = Number(fd.get("seed") || 0);
          const target = Number(fd.get("target") || -4);
          onStart(goal, seed, target);
        }}
      >
        <label className="block text-xs text-muted">
          Goal
          <input
            name="goal"
            required
            defaultValue="Find a strong binder in simulation mode"
            className="mt-1 w-full rounded border border-panel-border bg-white px-2 py-1.5 font-mono text-sm text-ink"
          />
        </label>
        <label className="block text-xs text-muted">
          Seed
          <input
            name="seed"
            type="number"
            defaultValue={0}
            className="mt-1 w-full rounded border border-panel-border bg-white px-2 py-1.5 font-mono text-sm"
          />
        </label>
        <label className="block text-xs text-muted">
          Target ≤
          <input
            name="target"
            type="number"
            step="0.1"
            defaultValue={-4}
            className="mt-1 w-full rounded border border-panel-border bg-white px-2 py-1.5 font-mono text-sm"
          />
        </label>
        <div className="flex items-end">
          <button
            type="submit"
            disabled={busy}
            className="w-full rounded bg-accent px-3 py-2 font-mono text-sm text-white disabled:opacity-50"
          >
            {busy ? "Running…" : "Start"}
          </button>
        </div>
      </form>
      {error ? (
        <p className="px-3 pb-2 font-mono text-xs text-red-800" data-testid="campaign-error">
          {error}
        </p>
      ) : null}
      {campaign ? (
        <div className="border-t border-panel-border px-3 py-2">
          <div className="mb-2 grid grid-cols-2 gap-2 font-mono text-[11px] md:grid-cols-4">
            <Stat label="status" value={campaign.status} />
            <Stat label="oracle calls" value={String(campaign.oracle_calls)} />
            <Stat
              label="best"
              value={campaign.best_value != null ? campaign.best_value.toFixed(3) : "—"}
            />
            <Stat label="labeled" value={String(campaign.n_labeled)} />
          </div>
          <table className="w-full font-mono text-[11px]" data-testid="iteration-table">
            <thead>
              <tr className="text-left text-muted">
                <th className="py-1">iter</th>
                <th>calls</th>
                <th>cost</th>
                <th>best</th>
                <th>acq</th>
              </tr>
            </thead>
            <tbody>
              {campaign.iterations.map((it) => (
                <tr key={it.iteration} className="border-t border-panel-border/70">
                  <td className="py-1">{it.iteration}</td>
                  <td>{it.oracle_calls}</td>
                  <td>{it.total_cost.toFixed(1)}</td>
                  <td>
                    {it.best_oracle_value != null ? it.best_oracle_value.toFixed(3) : "—"}
                  </td>
                  <td>{it.acquisition_method}</td>
                </tr>
              ))}
            </tbody>
          </table>
          <div className="mt-2">
            <ProvenanceBadge provenance={campaign.provenance} />
          </div>
        </div>
      ) : null}
    </section>
  );
}

function Stat({ label, value }: { label: string; value: string }) {
  return (
    <div className="rounded border border-panel-border bg-panel-inset px-2 py-1">
      <div className="text-muted">{label}</div>
      <div>{value}</div>
    </div>
  );
}
