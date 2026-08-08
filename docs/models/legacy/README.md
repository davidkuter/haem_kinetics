# Legacy model ladder (archived)

These pages and the matching modules under `haem_kinetics/models/legacy/` are the **previous** incremental ladder, kept for reference while the active ladder is rebuilt stepwise from Model 1.

| Legacy step | Change |
|-------------|--------|
| [Model 2](model2.md) | Model 1 + lipid factor φ on Fe(III) rates |
| [Model 3](model3.md) | `f_exp` uptake; enzymes track `f_exp`; φ on Hz |
| [Model 4](model4.md) | Aqueous ⇄ lipid Fe(III); Hz from lipid at full `k_hz` |
| [Model 5](model5.md) | Logistic enzyme clock (decoupled from `f_exp`) |
| [Model 6](model6.md) | Crystal-competent Fe(III) pool |

**Active ladder:** [docs/models.md](../models.md)

Import example:

```python
from haem_kinetics.models.legacy import LegacyModel3
```
