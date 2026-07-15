import { CampaignWorkspace } from "@/components/CampaignWorkspace";

export default function CampaignPage() {
  return (
    <div>
      <header className="border-b border-panel-border bg-white/80 px-4 py-2 backdrop-blur">
        <div className="mx-auto flex max-w-[1400px] items-baseline justify-between">
          <div className="flex items-baseline gap-3">
            <span className="font-mono text-xs uppercase tracking-[0.18em] text-accent">
              PeptideForge
            </span>
            <h1 className="font-display text-lg font-semibold text-ink">Campaign workspace</h1>
          </div>
          <span className="font-mono text-[10px] text-muted">
            simulation · provenance required · oracle gate open
          </span>
        </div>
      </header>
      <CampaignWorkspace />
    </div>
  );
}
