# Fractal Electrode Geometry & Scaling Laboratory

[![CI](https://github.com/AhmadEblagh205/sodium-ion-fractal-electrode/actions/workflows/ci.yml/badge.svg)](https://github.com/AhmadEblagh205/sodium-ion-fractal-electrode/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

This is a Python rebuild of my individual geometry and mathematics contribution to a UCL group design project on root-network-inspired sodium-ion electrode architectures. The wider project and its electrochemical design were group work; this repository covers my own section. The model is not a validated sodium-ion battery model, and the branch and electrochemical parameters are illustrative.

![Projected electrode render](reports/figures/electrode_render.png)

## What it does

- Generates a branching 2-D projected network with separate centreline shrink (lambda) and radius shrink (gamma).
- Measures the finite-depth box-counting slope and compares it with the ideal similarity expression log(b)/log(1/lambda).
- Computes surface and volume multipliers per generation and the radius window where surface grows while material falls.
- Implements De-Gennes rough-surface charge scaling and a branch-local charging-time relation.
- Exports a 2-D material field that conserves the cylindrical volume proxy.

## Results

For b=3 and lambda=0.65 the ideal similarity expression is about 2.55. The projected box-counting slope is not stable with depth: across ten random seeds its median rises from about 1.07 at depth 4 to 1.63 at depth 8. I report it as a finite-scale slope.

With b=3, lambda=0.65 and gamma=0.65, one ideal generation multiplies lateral surface by 1.2675 and material volume by 0.8239. Surface grows while volume falls for 0.5128 < gamma < 0.7161. A generation-5 branch has a local charging time of about 0.116 of the trunk value, which is a local scaling and not a whole-cell prediction.

At depth 7 the network has 1,093 segments, with a centreline length of about 103.3 trunk lengths. The 200x200 material field matches the analytic material proxy (0.04530123) to machine precision.

## How to run

```bash
pip install -r requirements.txt
PYTHONPATH=src python src/analyze.py
PYTHONPATH=src python src/export_field.py --depth 7 --resolution 200
pytest -q
```

The suite has 16 tests. Figures are written to `reports/figures/`.

## Notes

The model has no electrolyte transport, kinetics, thermal, mechanics or degradation component, so it says nothing about capacity or cycle life. My design idea was informed by Turon Teixidor, Park, Mukherjee, Kang and Madou, "Modeling fractal electrodes for Li-ion batteries", Electrochimica Acta 54 (2009) 5928-5936. That paper uses three-dimensional structures, and I use a simpler 2-D projected network. Derivations and the evidence boundary are in `docs/METHODOLOGY.md`.

## Licence

MIT, see [LICENSE](LICENSE).
