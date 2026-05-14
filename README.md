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

- [ ] data_input
- [ ] parse sample_saa_in_*.txt
- [ ] basic CLI execution

---

## NUMERICS

- [ ] set_x_array
- [ ] set_x_array_1d
- [ ] set_em_step_sizes
- [ ] dist_rebin

---

## PHYSICS

- [ ] photon_syn
- [ ] photon_ssc
- [ ] photon_alph
- [ ] calc_syn_fx
- [ ] calc_ssa_x_f
- [ ] calc_epse

---

## SOLVERS

- [ ] brents_minimum
- [x] newtons_method
- [ ] bw_euler_stepper

---

## VALIDATION

- [ ] Rust vs Fortran parity tests
- [ ] benchmark runtime
- [ ] benchmark memory usage
- [ ] verify interpolation equivalence