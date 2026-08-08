# haem_kinetics

Simulates the kinetics of haem speciation and haemozoin formation in the malaria
parasite *Plasmodium falciparum*.

## Models

Overview and comparison: **[docs/models.md](docs/models.md)**

Detailed pages (with process schematics):

| Model | Doc |
|-------|-----|
| Degradation | [docs/models/degradation.md](docs/models/degradation.md) |
| Model 1 | [docs/models/model1.md](docs/models/model1.md) |
| Model 2 | [docs/models/model2.md](docs/models/model2.md) |
| Model 3 | [docs/models/model3.md](docs/models/model3.md) |
| Model 4 | [docs/models/model4.md](docs/models/model4.md) |
| Model 5 | [docs/models/model5.md](docs/models/model5.md) |
| Model 6 | [docs/models/model6.md](docs/models/model6.md) |

## Quick start

```bash
pip install -e .
python examples/run.py
```

Edit `examples/run.py` to switch models and initial conditions.
