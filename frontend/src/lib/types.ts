/** Shared API types for the PeptideForge platform (P11). */

export type Provenance = {
  git_sha: string | null;
  data_version: string | null;
  tool_versions: Record<string, string>;
};

export type Iteration = {
  iteration: number;
  oracle_calls: number;
  total_cost: number;
  best_oracle_value: number | null;
  acquisition_method: string;
  status: string;
  notes: string | null;
};

export type Campaign = {
  campaign_id: string;
  status: string;
  goal: string;
  simulation_mode: boolean;
  reached_target: boolean;
  best_value: number | null;
  oracle_calls: number;
  total_cost: number;
  n_labeled: number;
  iterations: Iteration[];
  provenance: Provenance;
  agent_summary: string | null;
};

export type ParetoPoint = {
  candidate_id: string;
  sequence: string;
  neg_binding: number;
  solubility: number;
  oracle_value: number;
};

export type Structure = {
  candidate_id: string;
  sequence: string;
  pdb_text: string;
  fold_method: string;
  confidence: number;
  target_id: string;
  provenance: Provenance;
};

export type CalibrationBin = {
  bin_index: number;
  n: number;
  predicted_coverage: number;
  empirical_coverage: number;
  mean_interval_width: number;
};

export type Calibration = {
  n: number;
  coverage_target: number;
  empirical_coverage: number;
  ece: number;
  ece_threshold: number;
  passed: boolean;
  reliability_bins: CalibrationBin[];
  notes: string | null;
  provenance: Provenance;
};

export type TraceEvent = {
  kind: string;
  content: string;
  tool_name: string | null;
  data: Record<string, unknown>;
};

export type Trace = {
  session_id: string;
  events: TraceEvent[];
  summary: string | null;
  provenance: Provenance;
};

export type StartCampaignBody = {
  goal: string;
  seed: number;
  target_value: number;
  n_init: number;
  max_iterations: number;
  batch_size: number;
  simulation_mode: boolean;
  run_agent: boolean;
};
