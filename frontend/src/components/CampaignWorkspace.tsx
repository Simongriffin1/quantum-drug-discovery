/** Campaign workspace — DBTL iterations, Pareto, Mol*, calibration, agent trace. */

"use client";

import { useCallback, useState } from "react";
import { api, ApiError } from "@/lib/api";
import type {
  Calibration,
  Campaign,
  ParetoPoint,
  Structure,
  Trace,
} from "@/lib/types";
import { CalibrationPlot } from "@/components/CalibrationPlot";
import { CampaignForm } from "@/components/CampaignForm";
import { MolstarViewer } from "@/components/MolstarViewer";
import { ParetoPlot } from "@/components/ParetoPlot";
import { ReasoningTrace } from "@/components/ReasoningTrace";

export function CampaignWorkspace() {
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [campaign, setCampaign] = useState<Campaign | null>(null);
  const [pareto, setPareto] = useState<ParetoPoint[]>([]);
  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [structure, setStructure] = useState<Structure | null>(null);
  const [calibration, setCalibration] = useState<Calibration | null>(null);
  const [trace, setTrace] = useState<Trace | null>(null);

  const loadStructure = useCallback(async (campaignId: string, candidateId: string) => {
    const s = await api.getStructure(campaignId, candidateId);
    setStructure(s);
    setSelectedId(candidateId);
  }, []);

  const onStart = useCallback(
    async (goal: string, seed: number, target: number) => {
      setBusy(true);
      setError(null);
      try {
        const c = await api.startCampaign({
          goal,
          seed,
          target_value: target,
          n_init: 8,
          max_iterations: 3,
          batch_size: 2,
          simulation_mode: true,
          run_agent: true,
        });
        setCampaign(c);
        const [p, cal, tr] = await Promise.all([
          api.getPareto(c.campaign_id),
          api.getCalibration(c.campaign_id),
          api.getTrace(c.campaign_id),
        ]);
        setPareto(p);
        setCalibration(cal);
        setTrace(tr);
        if (p[0]) {
          await loadStructure(c.campaign_id, p[0].candidate_id);
        }
      } catch (e) {
        const msg = e instanceof ApiError ? e.message : "campaign failed";
        setError(msg);
      } finally {
        setBusy(false);
      }
    },
    [loadStructure],
  );

  const onSelect = useCallback(
    (id: string) => {
      if (!campaign) return;
      void loadStructure(campaign.campaign_id, id);
    },
    [campaign, loadStructure],
  );

  return (
    <div className="mx-auto max-w-[1400px] space-y-3 px-4 py-4" data-testid="campaign-workspace">
      <CampaignForm onStart={onStart} busy={busy} error={error} campaign={campaign} />
      <div className="grid gap-3 lg:grid-cols-2">
        <ParetoPlot
          points={pareto}
          selectedId={selectedId}
          onSelect={onSelect}
          provenance={
            campaign?.provenance ?? {
              git_sha: null,
              data_version: null,
              tool_versions: {},
            }
          }
        />
        <MolstarViewer structure={structure} />
        <CalibrationPlot calibration={calibration} />
        <ReasoningTrace trace={trace} />
      </div>
    </div>
  );
}
