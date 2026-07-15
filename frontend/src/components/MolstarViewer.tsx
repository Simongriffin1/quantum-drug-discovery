/** Mol* binding-mode viewer — loads molstar when available. */

"use client";

import { useEffect, useRef, useState } from "react";
import type { Structure } from "@/lib/types";
import { ProvenanceBadge } from "./ProvenanceBadge";

type Props = {
  structure: Structure | null;
};

export function MolstarViewer({ structure }: Props) {
  const hostRef = useRef<HTMLDivElement>(null);
  const [status, setStatus] = useState<string>("idle");

  useEffect(() => {
    const host = hostRef.current;
    if (!host || !structure?.pdb_text) {
      setStatus("empty");
      return;
    }
    let disposed = false;

    void (async () => {
      setStatus("loading");
      try {
        // Dynamic import keeps SSR clean; fail loud to status if missing.
        const molstar = await import("molstar/lib/viewer/viewer");
        if (disposed) return;
        host.innerHTML = "";
        const viewer = await molstar.Viewer.create(host, {
          layoutIsExpanded: false,
          layoutShowControls: false,
          layoutShowRemoteState: false,
          layoutShowSequence: true,
          layoutShowLog: false,
          layoutShowLeftPanel: false,
          viewportShowExpand: false,
          collapseLeftPanel: true,
        });
        await viewer.loadStructureFromData(structure.pdb_text, "pdb");
        if (!disposed) setStatus("ready");
      } catch (err) {
        // Lightweight CA ribbon fallback when molstar is not installed
        if (disposed || !host) return;
        host.innerHTML = "";
        const canvas = document.createElement("canvas");
        canvas.width = 480;
        canvas.height = 280;
        canvas.className = "mx-auto block bg-[#0c1216]";
        host.appendChild(canvas);
        drawCaFallback(canvas, structure.pdb_text);
        setStatus(
          `fallback (${err instanceof Error ? err.message.slice(0, 80) : "molstar unavailable"})`,
        );
      }
    })();

    return () => {
      disposed = true;
      if (host) host.innerHTML = "";
    };
  }, [structure]);

  return (
    <section className="panel" data-testid="molstar-viewer">
      <header className="panel-header">
        <h2>Binding mode (Mol*)</h2>
        <span className="panel-meta">{structure ? structure.fold_method : "—"}</span>
      </header>
      {!structure ? (
        <p className="px-3 py-8 text-sm text-muted">Select a Pareto point to load PDB.</p>
      ) : (
        <>
          <div className="px-3 py-1 font-mono text-[11px] text-muted">
            {structure.sequence} · conf={structure.confidence.toFixed(2)} · status={status}
          </div>
          <div ref={hostRef} className="h-[280px] w-full bg-[#0c1216]" data-testid="molstar-host" />
          <div className="panel-footer">
            <ProvenanceBadge provenance={structure.provenance} />
          </div>
        </>
      )}
    </section>
  );
}

/** Parse CA atoms and draw a simple connected backbone (not fabricated coords). */
export function drawCaFallback(canvas: HTMLCanvasElement, pdbText: string): void {
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  const atoms: { x: number; y: number; z: number }[] = [];
  for (const line of pdbText.split("\n")) {
    if (!line.startsWith("ATOM")) continue;
    if (line.slice(12, 16).trim() !== "CA") continue;
    const x = Number(line.slice(30, 38));
    const y = Number(line.slice(38, 46));
    const z = Number(line.slice(46, 54));
    if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) {
      atoms.push({ x, y, z });
    }
  }
  ctx.fillStyle = "#0c1216";
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  if (atoms.length === 0) {
    ctx.fillStyle = "#8fa3ad";
    ctx.font = "12px monospace";
    ctx.fillText("No CA atoms in PDB (synthetic stub ok)", 16, 24);
    return;
  }
  const xs = atoms.map((a) => a.x);
  const ys = atoms.map((a) => a.y);
  const minX = Math.min(...xs);
  const maxX = Math.max(...xs);
  const minY = Math.min(...ys);
  const maxY = Math.max(...ys);
  const dx = maxX - minX || 1;
  const dy = maxY - minY || 1;
  const pad = 24;
  const project = (a: { x: number; y: number }) => ({
    x: pad + ((a.x - minX) / dx) * (canvas.width - 2 * pad),
    y: pad + ((a.y - minY) / dy) * (canvas.height - 2 * pad),
  });
  ctx.strokeStyle = "#3d9b7a";
  ctx.lineWidth = 2;
  ctx.beginPath();
  atoms.forEach((a, i) => {
    const p = project(a);
    if (i === 0) ctx.moveTo(p.x, p.y);
    else ctx.lineTo(p.x, p.y);
  });
  ctx.stroke();
  ctx.fillStyle = "#7fd4b0";
  for (const a of atoms) {
    const p = project(a);
    ctx.beginPath();
    ctx.arc(p.x, p.y, 3, 0, Math.PI * 2);
    ctx.fill();
  }
}
