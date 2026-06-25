# semi-analytic-afterglow-rs

A Rust reimplementation of the original "semi-analytic-afterglow" Fortran project, aiming to keep the physics the same while making the code cleaner, safer, and (hopefully) faster.

Take a look at the [original Fortran code](https://github.fit.edu/GammaRayBurstAfterglow/semi-analytic-afterglow.git).

## Development Workflow

### 1. Pull the latest changes from `main`

```bash
git checkout main
git fetch origin
git pull origin main
```

### 2. Create a new feature branch

```bash
git checkout -b my-feature-branch
```

### 3. Make your changes

Build and run tests as needed:

```bash
cargo check
cargo test
```

### 4. Commit your work

```bash
git add .
git commit -m "describe your changes"
```

### 5. Push your branch to GitHub

```bash
git push -u origin my-feature-branch
```

Open a Pull Request from `my-feature-branch` into `main`.

### 6. Sync changes to the FIT repository

This project maintains a synchronized mirror on the FIT GitHub Enterprise server.

After committing your changes, run:

```bash
./scripts/push-fit.sh
```

This script will:

1. Create/reset the local `fit-runners` branch from `fit/main`
2. Cherry-pick commits from your current branch
3. Push the result to the FIT synchronization branch
4. Return you to your original branch

### Updating an Existing Branch

If `main` has moved forward since your branch was created:

```bash
git checkout main
git pull origin main

git checkout my-feature-branch
git rebase main
```

Resolve any conflicts, then continue:

```bash
git rebase --continue
```

Push the updated branch:

```bash
git push --force-with-lease
```

---

## FOUNDATION

- [x] parameters.rs
- [x] constants.rs
- [x] inputs.rs
- [x] globals.rs
- [x] opacity.rs

---

## FUNCTIONS

- [x] E_of_theta
- [x] Gam_of_theta
- [x] R_of_theta
- [x] theta_f
- [x] d_rperp_dmu
- [x] fixed_rperp
- [x] gps_egg
- [x] r_theta_quartic

---

## INPUT / DRIVER

- [x] data_input
- [ ] parse sample_saa_in_*.txt
- [ ] basic CLI execution

---

## NUMERICS

- [x] set_x_array
- [x] set_x_array_1d
- [x] set_x_array_simpsons
- [x] set_em_step_sizes
- [x] dist_rebin

---

## PHYSICS

- [x] photon_syn
- [x] photon_ssc
- [x] photon_alph
- [x] calc_syn_fx
- [x] calc_ssa_x_f
- [x] calc_epse
- [x] set_elec_dist
- [x] rad_transfer_diffeq

---

## SOLVERS / SPECIAL FUNCTIONS

- [x] brents_minimum
- [x] newtons_method
- [x] bw_euler_stepper
- [x] bessik
- [x] beschb
- [x] chebev_subr

---

## VALIDATION

- [ ] Rust vs Fortran parity tests
- [ ] benchmark runtime
- [ ] benchmark memory usage
- [ ] verify interpolation equivalence
