import { describe, expect, it } from "vitest";
import { drawCaFallback } from "@/components/MolstarViewer";
import type { Calibration, ParetoPoint, Provenance, Trace } from "@/lib/types";

const mockProv: Provenance = {
  git_sha: "abc12345deadbeef",
  data_version: "synthetic_v1",
  tool_versions: { peptideforge: "0.1.0", platform: "p11" },
};

const mockPareto: ParetoPoint[] = [
  {
    candidate_id: "c1",
    sequence: "AAAAA",
    neg_binding: 1.2,
    solubility: 0.5,
    oracle_value: -1.2,
  },
  {
    candidate_id: "c2",
    sequence: "FLIVV",
    neg_binding: 2.0,
    solubility: 0.1,
    oracle_value: -2.0,
  },
];

const mockCal: Calibration = {
  n: 8,
  coverage_target: 0.9,
  empirical_coverage: 0.875,
  ece: 0.025,
  ece_threshold: 0.1,
  passed: true,
  reliability_bins: [
    {
      bin_index: 0,
      n: 8,
      predicted_coverage: 0.9,
      empirical_coverage: 0.875,
      mean_interval_width: 1.0,
    },
  ],
  notes: "mock",
  provenance: mockProv,
};

const mockTrace: Trace = {
  session_id: "sess-1",
  events: [
    { kind: "gate_pause", content: "oracle_validity skipped", tool_name: null, data: {} },
    { kind: "tool_result", content: "campaign done", tool_name: "run_simulation_campaign", data: {} },
  ],
  summary: "best_oracle_value = -2.0",
  provenance: mockProv,
};

describe("platform mock artifacts", () => {
  it("pareto points are selectable by id", () => {
    expect(mockPareto.map((p) => p.candidate_id)).toContain("c1");
  });

  it("calibration gate fields are present", () => {
    expect(mockCal.ece).toBeLessThan(mockCal.ece_threshold);
    expect(mockCal.provenance.data_version).toBe("synthetic_v1");
  });

  it("trace includes gate_pause and tool_result", () => {
    const kinds = mockTrace.events.map((e) => e.kind);
    expect(kinds).toContain("gate_pause");
    expect(kinds).toContain("tool_result");
  });

  it("drawCaFallback renders without throwing", () => {
    const canvas = {
      width: 100,
      height: 80,
      getContext: () => ({
        fillStyle: "",
        strokeStyle: "",
        lineWidth: 0,
        font: "",
        fillRect: () => undefined,
        fillText: () => undefined,
        beginPath: () => undefined,
        moveTo: () => undefined,
        lineTo: () => undefined,
        stroke: () => undefined,
        arc: () => undefined,
        fill: () => undefined,
      }),
    } as unknown as HTMLCanvasElement;
    const pdb =
      "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 50.00           C\n" +
      "ATOM      2  CA  ALA A   2       1.000   1.000   0.000  1.00 50.00           C\n";
    expect(() => drawCaFallback(canvas, pdb)).not.toThrow();
  });
});
