"""ESM-2 masked language-model sampling for peptide proposals.

Method: randomly mask positions in a seed sequence, run ESM-2 forward pass,
sample amino acids at masked sites from the LM head (Rives et al. / Lin et al.,
ESM-2). Requires ``fair-esm`` + PyTorch — missing deps raise ``ESM2UnavailableError``.
"""

from __future__ import annotations

import random
from typing import Any

from peptideforge.contracts.models import Candidates, PeptideCandidate
from peptideforge.generators.deps import ESM2UnavailableError, require_esm2
from peptideforge.generators.filters import (
    is_diverse_enough,
    is_valid_sequence,
    normalize_allowed_residues,
)
from peptideforge.generators.mutation import _parse_constraints

# Small public checkpoint for optional local runs; override via ESM2Generator(model_name=...).
DEFAULT_MODEL = "esm2_t6_8M_UR50D"
MODEL_LOADERS: dict[str, str] = {
    "esm2_t6_8M_UR50D": "esm2_t6_8M_UR50D",
    "esm2_t12_35M_UR50D": "esm2_t12_35M_UR50D",
}


class ESM2Generator:
    """Masked LM sampling with Meta ESM-2 (fair-esm)."""

    def __init__(
        self,
        *,
        model_name: str = DEFAULT_MODEL,
        masks_per_sample: tuple[int, int] = (1, 4),
        max_identity: float = 0.95,
        generation_method: str = "esm2_masked",
        device: str | None = None,
    ) -> None:
        if model_name not in MODEL_LOADERS:
            raise ValueError(
                f"unsupported ESM-2 model {model_name!r}; choose from {sorted(MODEL_LOADERS)}"
            )
        self.model_name = model_name
        lo, hi = masks_per_sample
        if lo < 1 or hi < lo:
            raise ValueError("masks_per_sample must be (lo, hi) with 1 <= lo <= hi")
        self.masks_per_sample = masks_per_sample
        self.max_identity = max_identity
        self.generation_method = generation_method
        self.device = device
        self._bundle: tuple[Any, Any, Any, Any] | None = None

    def propose(
        self,
        *,
        n: int,
        seed_sequences: tuple[str, ...] | None = None,
        seed: int | None = None,
        constraints: dict[str, object] | None = None,
    ) -> Candidates:
        if n < 1:
            raise ValueError("n must be >= 1")
        if seed_sequences is None or not seed_sequences:
            raise ValueError("ESM2Generator requires seed_sequences")
        rng = random.Random(seed)
        cfg = _parse_constraints(constraints, default_allowed=normalize_allowed_residues(None))
        model, alphabet, batch_converter, torch = self._load()
        seeds = [s.upper() for s in seed_sequences]
        for seq in seeds:
            if not is_valid_sequence(
                seq,
                min_length=cfg.min_length,
                max_length=cfg.max_length,
                allowed_residues=cfg.allowed_residues,
                fixed_positions=cfg.fixed_positions_dict,
            ):
                raise ValueError(f"invalid seed sequence: {seq!r}")

        selected: list[str] = []
        items: list[PeptideCandidate] = []
        attempts = 0
        max_attempts = n * 300

        while len(items) < n and attempts < max_attempts:
            attempts += 1
            parent = rng.choice(seeds)
            child = self._sample_masked(
                parent,
                rng,
                cfg,
                model=model,
                alphabet=alphabet,
                batch_converter=batch_converter,
                torch=torch,
            )
            if child is None:
                continue
            if not is_diverse_enough(child, selected, max_identity=self.max_identity):
                continue
            selected.append(child)
            items.append(
                PeptideCandidate(
                    sequence=child,
                    generation_method=self.generation_method,
                    metadata={"parent_sequence": parent, "esm2_model": self.model_name},
                )
            )

        if len(items) < n:
            raise RuntimeError(
                f"ESM2Generator produced {len(items)}/{n} diverse valid sequences "
                f"after {attempts} attempts"
            )
        return Candidates(items=tuple(items), seed=seed)

    def _load(self) -> tuple[Any, Any, Any, Any]:
        if self._bundle is not None:
            return self._bundle
        esm, torch = require_esm2()
        loader_name = MODEL_LOADERS[self.model_name]
        loader = getattr(esm.pretrained, loader_name)
        model, alphabet = loader()
        device = self.device or ("cuda" if torch.cuda.is_available() else "cpu")
        model = model.to(device)
        model.eval()
        batch_converter = alphabet.get_batch_converter()
        self._bundle = (model, alphabet, batch_converter, torch)
        return self._bundle

    def _sample_masked(
        self,
        parent: str,
        rng: random.Random,
        cfg: Any,
        *,
        model: Any,
        alphabet: Any,
        batch_converter: Any,
        torch: Any,
    ) -> str | None:
        seq = list(parent)
        n_masks = rng.randint(*self.masks_per_sample)
        mutable = [
            i
            for i in range(len(seq))
            if i not in cfg.fixed_positions and seq[i] in cfg.allowed_residues
        ]
        if not mutable:
            return None
        mask_positions = rng.sample(mutable, min(n_masks, len(mutable)))
        masked_tokens = list(seq)
        for pos in mask_positions:
            masked_tokens[pos] = "<mask>"
        _, _, batch_tokens = batch_converter([("peptide", "".join(masked_tokens))])
        device = next(model.parameters()).device
        batch_tokens = batch_tokens.to(device)

        mask_idx = alphabet.mask_idx
        with torch.no_grad():
            out = model(batch_tokens, repr_layers=[], return_contacts=False)
            logits = out["logits"][0]

        for pos in mask_positions:
            # +1 for BOS token in ESM batch format
            token_logits = logits[pos + 1]
            probs = torch.softmax(token_logits, dim=-1)
            # Sample until allowed residue chosen
            for _ in range(50):
                sampled_id = int(torch.multinomial(probs, 1).item())
                aa = alphabet.get_tok(sampled_id)
                if len(aa) == 1 and aa in cfg.allowed_residues and aa != seq[pos]:
                    seq[pos] = aa
                    break
            else:
                return None

        candidate = "".join(seq)
        if not is_valid_sequence(
            candidate,
            min_length=cfg.min_length,
            max_length=cfg.max_length,
            allowed_residues=cfg.allowed_residues,
            fixed_positions=cfg.fixed_positions_dict,
        ):
            return None
        return candidate


def check_esm2_available() -> None:
    """Raise ``ESM2UnavailableError`` if fair-esm is not importable."""
    require_esm2()
