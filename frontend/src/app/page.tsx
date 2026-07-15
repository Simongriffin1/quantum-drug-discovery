import Link from "next/link";

export default function HomePage() {
  return (
    <main className="mx-auto flex min-h-[70vh] max-w-3xl flex-col justify-center gap-6 px-6 py-16">
      <p className="font-mono text-xs uppercase tracking-[0.2em] text-accent">PeptideForge</p>
      <h1 className="font-display text-4xl font-semibold tracking-tight text-ink md:text-5xl">
        Physics-grounded peptide design
      </h1>
      <p className="max-w-xl text-lg leading-relaxed text-muted">
        Closed-loop campaigns driven by an MM-GBSA / quantum chemistry oracle — no proprietary
        training data. Open the workspace to run a synthetic-mode campaign and inspect Pareto,
        Mol*, calibration, and agent trace with full provenance.
      </p>
      <Link
        href="/campaign"
        className="inline-flex w-fit rounded bg-accent px-4 py-2 font-mono text-sm text-white"
      >
        Open campaign workspace
      </Link>
    </main>
  );
}
