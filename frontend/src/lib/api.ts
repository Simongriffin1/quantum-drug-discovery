/** Client for PeptideForge backend API. */

import type {
  Calibration,
  Campaign,
  ParetoPoint,
  StartCampaignBody,
  Structure,
  Trace,
} from "./types";

const DEFAULT_BASE = process.env.NEXT_PUBLIC_API_URL ?? "http://127.0.0.1:8000";

export class ApiError extends Error {
  constructor(
    message: string,
    readonly status: number,
  ) {
    super(message);
  }
}

async function request<T>(path: string, init?: RequestInit, base = DEFAULT_BASE): Promise<T> {
  const res = await fetch(`${base}${path}`, {
    ...init,
    headers: {
      "Content-Type": "application/json",
      ...(init?.headers ?? {}),
    },
  });
  if (!res.ok) {
    const detail = await res.text();
    throw new ApiError(detail || res.statusText, res.status);
  }
  return res.json() as Promise<T>;
}

export function createApi(baseUrl = DEFAULT_BASE) {
  return {
    startCampaign: (body: StartCampaignBody) =>
      request<Campaign>("/campaigns", { method: "POST", body: JSON.stringify(body) }, baseUrl),
    getCampaign: (id: string) => request<Campaign>(`/campaigns/${id}`, undefined, baseUrl),
    getPareto: (id: string) => request<ParetoPoint[]>(`/campaigns/${id}/pareto`, undefined, baseUrl),
    getStructure: (id: string, candidateId: string) =>
      request<Structure>(`/campaigns/${id}/structures/${candidateId}`, undefined, baseUrl),
    getCalibration: (id: string) =>
      request<Calibration>(`/campaigns/${id}/calibration`, undefined, baseUrl),
    getTrace: (id: string) => request<Trace>(`/campaigns/${id}/trace`, undefined, baseUrl),
  };
}

export const api = createApi();
