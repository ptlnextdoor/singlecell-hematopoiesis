# asiri-domain-template

Template anchor repo for the Round 3 (Stanford / Asiri) outreach. **Specialize
this once the target's research area is known.**

## How to specialize (do this before Round 3)

1. Identify Asiri's full name, lab, and subfield (e.g. NLP, systems, vision, HCI,
   theory). Aayushya must provide this.
2. Pick ONE recent paper from that lab.
3. Replace this README and `src/` with a small, honest reproduction or evaluation
   harness for a specific claim in that paper, following the pattern in the sibling
   repos:
   - a `scripts/selfcheck.py` that runs with no large download and asserts correctness
   - a real result (one plot or one metric table)
   - a stub adapter for the real public dataset
4. Keep it a capability repo reused across every email in the Round 3 cluster;
   put the paper-specific hook in each email, not in a new repo per person.

## Pattern to copy

See `../eeg-physio-ml` (a full working example: synthetic-data self-check,
frozen-feature probes, dataset stub) and `../search-planning-robotics`
(algorithm reproduction with a correctness self-check).

## Status

Placeholder. Blocked on target field.

## License

MIT
