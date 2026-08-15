# haem_kinetics

Simulates the kinetics of haem speciation and haemozoin formation in the malaria
parasite *Plasmodium falciparum*.

## Models

Overview and comparison: **[docs/models.md](docs/models.md)**

Garnie Hb/Hm/Hz assay vs pHrodo: **[docs/garnie_fractionation.md](docs/garnie_fractionation.md)**

Detailed pages (with process schematics):

| Model | Doc |
|-------|-----|
| Degradation | [docs/models/degradation.md](docs/models/degradation.md) |
| Model 1 | [docs/models/model1.md](docs/models/model1.md) |
| Model 2a / 2b | [docs/models/model2.md](docs/models/model2.md) |
| Model 3 | [docs/models/model3.md](docs/models/model3.md) |
| Model 4a / 4b | [docs/models/model4.md](docs/models/model4.md) |
| Model 5 | [docs/models/model5.md](docs/models/model5.md) |
| Model 6 | [docs/models/model6.md](docs/models/model6.md) |
| Model 7 | [docs/models/model7.md](docs/models/model7.md) |
| Model 8 | [docs/models/model8.md](docs/models/model8.md) |
| Model 9a / 9b / 9c | [docs/models/model9.md](docs/models/model9.md) |
| Model 10 | [docs/models/model10.md](docs/models/model10.md) |
| Model 99 (what-if) | [docs/models/model99.md](docs/models/model99.md) |
| Legacy 2–6 | [docs/models/legacy/](docs/models/legacy/README.md) |

## Quick start

```bash
pip install -e .
python examples/run.py
```

Edit `examples/run.py` to switch models and initial conditions.
