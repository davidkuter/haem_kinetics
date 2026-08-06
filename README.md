# haem_kinetics

Simulates the kinetics of haem speciation and haemozoin formation in the malaria
parasite *Plasmodium falciparum*.

## Models

See **[docs/models.md](docs/models.md)** for a full comparison of Models 1–7 and the
Degradation sandbox: assumptions, state variables, and how each version differs.

| Model | Role |
|-------|------|
| 1–4, Degradation | Legacy ladder (linear → exponential uptake, lipid `φ` on Hz) |
| **5** | Corrected aqueous/lipid Fe(III) pools; Hz from lipid at `k_hz` |
| **6** | + dynamic DV volume, depleting host Fe, sigmoidal uptake (fg/cell) |
| **7** | + falcipain-2/3 |

## Quick start

```bash
pip install -e .
python examples/run.py
```

Edit `examples/run.py` to switch models and initial conditions.
