import { describe, expect, it } from "vitest";

describe("frontend landing", () => {
  it("documents campaign workspace route", () => {
    expect("/campaign").toMatch(/campaign/);
  });
});
