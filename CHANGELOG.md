# Changelog

## 0.3.0, defensibility audit

- Separated the 2-D projected box-counting slope, ideal similarity dimension and 3-D electrochemical De-Gennes surface dimension instead of treating them as interchangeable.
- Retired the unsupported “optimum fractal dimension = 2.55” wording; 2.55 is now reported only as an ideal similarity scaling quantity for `b=3`, `lambda=0.65`.
- Separated centreline-length shrink from branch-radius shrink and replaced plotting linewidth as a pseudo-physical width with an explicit cylindrical radius proxy.
- Derived and tested the per-generation cylindrical surface and material-volume multipliers `b*lambda*gamma` and `b*lambda*gamma^2`.
- Added the radius-shrink window where ideal lateral surface grows while added material volume falls; the default 0.65 lies inside the interval 0.5128–0.7161.
- Added a ten-seed finite-depth box-counting stability study; the projected slope rises strongly with branching depth and is no longer presented as a depth-invariant material constant.
- Replaced generation-number “transport depth” with a material-volume-weighted connected root-path metric.
- Tightened domains and validation for the De-Gennes and branch-local charging-time relations.
- Rebuilt the density exporter around cylindrical material volume with conservative bilinear deposition plus centroid-fidelity diagnostics.
- Removed the claim that the exported 2-D field is essentially ready to plug directly into PyBaMM.
- Added `docs/METHODOLOGY.md` and expanded tests from 6 to 16.

## 0.2.0, update round 3

Added a density field exporter: the grown network rasterised to a material-per-cell grid, saved as npz with coordinates plus a rendered field, the bridge toward a heterogeneous electrode simulation.

## 0.1.0

Initial release.
