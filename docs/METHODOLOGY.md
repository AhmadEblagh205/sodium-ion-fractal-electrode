# Methodology and evidence boundary

## Scope

This repository is a clean Python rebuild of my individual geometry/mathematical contribution from a UCL group design project on root-network-inspired sodium-ion electrode architectures. It is a **design-study / scientific-computing project**, not a validated sodium-ion cell model.

The audit deliberately separates three ideas that the earlier version blurred together:

1. a **2-D projected branching skeleton** used for geometry and visualisation;
2. **3-D cylindrical branch proxies** used to reason about lateral surface area and material volume;
3. **electrochemical fractal-surface scaling** such as the De-Gennes relation, where the rough-surface dimension is conventionally between 2 and 3.

Those dimensions are not interchangeable.

## Geometry

A branch has centreline length `L` and cylindrical radius `r`. A parent creates `b` children. Centreline length shrinks by `lambda` and radius by `gamma` at each generation. The visual layout is a stochastic 2-D projection, while the surface/volume proxies are calculated from the cylindrical branch dimensions.

For an ideal deterministic generation, ignoring end caps and overlaps,

- lateral surface scales as `b * lambda * gamma`;
- material volume scales as `b * lambda * gamma^2`.

Therefore a useful **surface-efficient branching window** exists when

`b*lambda*gamma > 1` and `b*lambda*gamma^2 < 1`.

For the default `b=3`, `lambda=0.65`, this becomes

`0.5128 < gamma < 0.7161`.

The default `gamma=0.65` lies inside that interval, giving ideal per-generation multipliers of 1.2675 for lateral area and 0.823875 for added material volume. This is a scaling result, not a battery-performance optimum.

## Similarity dimension vs projected box counting

`log(b)/log(1/lambda)` is reported only as an **ideal similarity dimension**. For `b=3`, `lambda=0.65` it is about 2.55. A non-overlapping self-similar set embedded in a plane cannot realise a dimension above 2, so that number must not be presented as the fractal dimension of the rendered 2-D tree.

The box-counting calculation is therefore reported as a **finite-scale projected box-counting slope**. It depends strongly on finite branching depth and sampled scale range. In the checked-in 10-seed study the median slope rises from about 1.07 at depth 4 to about 1.63 at depth 8. This is evidence that a single finite-tree slope should not be advertised as a depth-invariant material constant.

## Electrochemical scaling

Two preliminary relations from the design-study context are retained as illustrative scalings.

### De-Gennes rough-surface charge scaling

`Q(t) ~ t^((3-D)/2)` for a rough surface dimension `2 <= D <= 3`.

The code explicitly validates that domain. The `D` in this relation describes a rough electrochemical surface in three-dimensional space. It is not the projected box slope of the 2-D branching drawing.

### Branch-local reaction-zone charging time

The branch-axis time scale is

`t_i = r0 * gamma^i * (1-eps) * q_plus / (2*j)`.

Hence `t_i/t_0 = gamma^i`. For `gamma=0.65`, generation 5 has a local time scale of about 0.116 of the trunk value under this model. This does **not** predict whole-electrode charge time; network coupling, transport, electrochemistry and cell-scale constraints are omitted.

## Conservative material field

The exporter now rasterises the cylindrical material-volume proxy `pi*r^2*L`, not a plotting linewidth. Material is deposited bilinearly on a regular 2-D grid. Two invariants are checked:

- the field integral equals the analytical cylindrical-volume proxy to floating-point precision;
- the raster centre of mass matches the analytical segment-weighted centre of mass to within a grid-cell scale.

The resulting field is a generic spatial representation for downstream modelling. It is **not** described as a drop-in PyBaMM input; a separate model mapping would be required.

## Main limitations

- The geometry is a 2-D projection of a concept motivated by 3-D branching electrodes.
- Cylindrical branch surface/volume quantities are geometric proxies.
- Electrochemical parameters are illustrative, not fitted to sodium-ion experiments.
- There is no electrolyte transport, solid diffusion, potential field, kinetics, thermal model, mechanical constraint or degradation model.
- No claim is made that the default branching parameters are optimal for a real battery.

## Literature context

The original design idea was informed by work such as G. Turon Teixidor *et al.*, “Modeling fractal electrodes for Li-ion batteries,” *Electrochimica Acta* 54 (2009) 5928–5936, DOI: 10.1016/j.electacta.2009.05.060. That work models three-dimensional non-deterministic branching structures composed of cylindrical elements. The current repository uses a simpler 2-D projected geometry and is explicit about that difference.

## Notes moved from the README

Three quantities had been conflated in the original repository: the box-counting slope of a finite 2-D projected tree, the ideal self-similar quantity log(b)/log(1/lambda), and the De-Gennes dimension D in Q ~ t^((3-D)/2), which is conventionally defined for a surface embedded in 3-D. For b=3 and lambda=0.65 the ideal expression is about 2.55. That can be a 3-D scaling quantity but cannot be the dimension of a non-overlapping self-similar set in the rendered plane. Across ten random seeds the median projected box-counting slope rises from about 1.07 at depth 4 to 1.63 at depth 8, so I report it as a finite-scale slope and not a material constant.

Default surface-efficient branching (b=3, lambda=0.65, gamma=0.65): surface multiplier 1.2675 and volume multiplier 0.8239 per generation, that is 26.75% more lateral area and 17.61% less material volume. This is a geometric result and does not claim gamma=0.65 is experimentally optimal.

Finite network at depth 7: 1,093 segments, centreline length about 103.3 trunk lengths, material-volume-weighted connected path about 1.60 trunk lengths. The old code used a visual line width as a physical width, and v0.3 replaced it with a dimensionless cylindrical radius proxy.

Material field. The default 200x200 export gives an analytic material proxy of 0.04530123 and an identical field integral, so the relative conservation error is at machine precision. The raster centre of mass is also checked against the analytical segment-weighted centre of mass. An earlier claim that the field was ready for PyBaMM was removed, since mapping a 2-D architecture into a battery PDE model is a separate task.

Not established: battery capacity, rate capability, cycle life, an experimentally optimal fractal dimension or a validated sodium-ion electrode design. There is no transport PDE, kinetics, thermal, mechanics or degradation model.

Test suite (16 tests): box counting on a straight line, branch-count and geometry invariants, the 2.55 similarity calculation and its planar caveat, surface and volume multipliers, the radius window, geometry versus analytic agreement, finite-depth box-slope sensitivity, De-Gennes limits and domain validation, branch-local charging time, field conservation, and raster-centroid fidelity and resolution invariance.
