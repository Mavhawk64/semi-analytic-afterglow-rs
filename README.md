# semi-analytic-afterglow-rs

A Rust reimplementation of the original "semi-analytic-afterglow" Fortran project, aiming to keep the physics the same while making the code cleaner, safer, and (hopefully) faster.

Take a look at the [original Fortran code](https://github.fit.edu/GammaRayBurstAfterglow/semi-analytic-afterglow.git).

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
