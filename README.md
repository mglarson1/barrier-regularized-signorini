# Barrier-Regularized Signorini Reproducibility Code

MATLAB implementation and reproducibility package for the numerical examples in
*A Barrier-Regularized Symmetric Nitsche Method for the Signorini Problem* by
Peter Hansbo and Mats G. Larson.

The code solves the scalar Signorini problem with linear and quadratic finite
elements, a symmetric Nitsche formulation, and a smooth logarithmic-barrier
regularization. It contains the experiment drivers, saved reference datasets,
and figure-generation scripts used for the paper.

## Requirements

- MATLAB (tested with MATLAB R2025a)
- No additional MATLAB toolboxes

MATLAB's dependency analysis reports only the base MATLAB product for the main
reproducibility script.

## Quick check

Start MATLAB, change to the repository directory, and run:

```matlab
check_barrier_evaluation
```

This verifies the stable, cancellation-free evaluation of the barrier term and
its derivative.

## Complete reproduction

Run:

```matlab
run_scalar_reproducibility
```

The script:

1. checks the stable barrier evaluation;
2. regenerates all saved numerical datasets;
3. writes the paper figures to `figures/`.

The complete run contains mesh-refinement studies and parameter sweeps and can
take substantially longer than the quick check.

## Experiments and outputs

- `run_reg_sweep.m` and `run_reg_sweep2.m` study the regularization error
  relative to the unregularized discrete solution. They write `results_p1.mat`
  and `results_p2.mat`.
- `run_exact_all.m` performs convergence studies for exact solutions A and B
  and writes `exact_conv_data.mat`.
- `run_numerics_extra.m` records Newton histories, contact-pressure error,
  penetration and feasibility diagnostics, and an artificial-gap sweep. It
  writes `numerics_extra.mat`.
- `run_numerics_gamma_mesh.m` performs the Nitsche-parameter and unstructured-
  mesh studies and writes `numerics_gamma_mesh.mat`.
- Scripts whose names begin with `make_` regenerate the figures from the saved
  datasets.

The tracked `.mat` files contain the reference results used to generate the
paper figures. Rerunning the experiments overwrites them, making it possible to
compare newly generated results with the repository versions using Git.

## Numerical conventions

The implementation sets `epsilon = dl^pow`, where `dl` is the contact-edge
length, and evaluates

```text
sqrt(w^2 + epsilon) = sqrt(w^2 + 4s),   s = epsilon/4.
```

Thus, on the uniform meshes, `s = h^pow/4`. The barrier value and derivative
use a cancellation-free branch for negative `w`.

The Nitsche parameter is local, `gamma = gamma0/dl`. The structured studies use
`gamma0 = 10` for P1 and `gamma0 = 20` for P2 unless the parameter is being
swept. The unstructured meshes retain contact-boundary edges of length `1/n`,
so the same relation `s = h^pow/4` holds on the contact boundary.

The regularization sweeps report the `H1` seminorm of the difference between
two discrete solutions. The exact-solution studies report the full `H1` norm.
P1 model-problem assembly uses one-point volume and two-point edge quadrature;
P2 model-problem assembly uses degree-five volume and four-point edge
quadrature. The exact-solution studies use degree-five volume and four-point
edge quadrature for both element orders.

The model-problem solvers stop when the Euclidean residual on the free degrees
of freedom is below `max(1e-10,1e-12*r0)`. The exact-solution solver uses
`max(1e-11,1e-12*r0)`. All reported runs use undamped Newton steps.

## Citation

If you use this software, please cite the accompanying paper. GitHub can also
generate software citation metadata from `CITATION.cff`. Publication details
and a DOI can be added when they become available.

## License

This project is available under the BSD 3-Clause License. See `LICENSE`.
