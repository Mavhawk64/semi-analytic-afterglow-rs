use crate::modules::globals::Globals;

/// Goes from x_start to x_end using the backwards (implicit) Euler method.
/// Because we know a priori that our function varies most rapidly near
/// x_end, we reduce the step size periodically to enhance accuracy while
/// conserving computation time.
///
/// # Arguments
///
/// * `func` - external function computing the value of f(x,y)
/// * `x_start` - initial value of x
/// * `x_end` - final value of x. Note that x_end < x_start is allowed
/// * `y_start` - initial values of y
/// * `r_params` - array of up to 7 floating point parameters to be used by func
/// * `i_params` - array of up to 7 integer parameters to be used by func
/// * `globals` - runtime global state containing `em_steps` and `em_size_array`
pub fn bw_euler_stepper<F>(
    func: F,
    x_start: f64,
    x_end: f64,
    y_start: &[f64],
    r_params: &[f64; 7],
    i_params: &[i32; 7],
    globals: &Globals,
) -> Vec<f64>
where
    F: Fn(f64, usize, &mut [f64], &mut [f64], &[f64; 7], &[i32; 7]),
{
    let n_y = y_start.len();

    let x_range = x_end - x_start;

    // current location
    let mut x_curr = x_start;

    let mut y_curr = y_start.to_vec();

    let mut i_params_int = *i_params;
    let mut r_params_int = *r_params;

    let mut f_out1 = vec![0.0; n_y];
    let mut f_out2 = vec![0.0; n_y];

    // The ftemp arrays are only used in the midpoint evaluation
    let mut ftemp1 = vec![0.0; n_y];
    let mut ftemp2 = vec![0.0; n_y];

    //--------------------------------------------------------------------------
    // For the first step only (getting y_1), use the midpoint method to cut
    // down on error.
    //
    // This requires evaluating func at both the start and the midpoint of
    // the first step. (need x_0 and x_1/2 -- never x_1)
    //--------------------------------------------------------------------------

    // Evaluate at the start of the step
    // (size of the current step)
    let mut del_x = globals.em_size_array[0] * x_range;

    i_params_int[2] = 1; // 1 = i_tm, 2 = k_x, and 3 = n_step
    i_params_int[3] = 10; // 0 -> s not equal to s_near or s_far; 10 -> s equal

    r_params_int[6] = 0.5 * del_x;

    func(
        x_curr,
        n_y,
        &mut ftemp1,
        &mut ftemp2,
        &r_params_int,
        &i_params_int,
    );

    // Evaluate at the midpoint
    i_params_int[3] = 0; // 0 -> s not equal to s_near or s_far; 10 -> s equal

    r_params_int[6] = del_x;

    func(
        x_curr + 0.5 * del_x,
        n_y,
        &mut f_out1,
        &mut f_out2,
        &r_params_int,
        &i_params_int,
    );

    // Combine everything into the evaluation of y_next
    for i in 0..n_y {
        y_curr[i] = y_curr[i] + f_out1[i] * del_x
            - f_out2[i] * del_x * (y_curr[i] + 0.5 * del_x * (ftemp1[i] - y_curr[i] * ftemp2[i]));
    }

    // And make sure we update x_curr before entering the loop
    x_curr += del_x;

    //--------------------------------------------------------------------------
    // Loop over steps in the backwards Euler method
    //--------------------------------------------------------------------------

    for n_step in 2..=globals.em_steps {
        // We reduce the logic from the original code
        // because it was redundant.
        del_x = globals.em_size_array[n_step - 1] * x_range;

        // This has caused some confusion.
        // The reason why we update x_curr here rather than at the end is
        // because we update then calculate, and y_n depends on x_n,
        // not x_{n-1}.
        // We start at y_2 here, since y_0 is given and y_1 is calculated
        // with the midpoint method above.
        x_curr += del_x;

        i_params_int[2] = n_step as i32; // 1 = i_tm, 2 = k_x, and 3 = n_step
        i_params_int[3] = 0; // 0 -> s not equal to s_near or s_far; 10 -> s equal

        r_params_int[6] = del_x;

        func(
            x_curr,
            n_y,
            &mut f_out1,
            &mut f_out2,
            &r_params_int,
            &i_params_int,
        );

        for i in 0..n_y {
            y_curr[i] = (y_curr[i] + f_out1[i] * del_x) / (1.0 + f_out2[i] * del_x);
        }
    }

    //--------------------------------------------------------------------------
    // For the last step, use the midpoint method rather than the backwards
    // Euler method.
    //
    // Evaluating solely at the final point leads to an unrealistic photon
    // production from the highest-energy electrons, which would be completely
    // uncooled since the final point lies on the shock itself.
    //--------------------------------------------------------------------------

    // Evaluate at the start of the step
    del_x = x_end - x_curr;

    i_params_int[2] = globals.em_steps as i32 + 1; // 1 = i_tm, 2 = k_x, and 3 = n_step

    i_params_int[3] = 0; // 0 -> s not equal to s_near or s_far; 10 -> s equal

    r_params_int[6] = 0.5 * del_x;

    func(
        x_curr,
        n_y,
        &mut ftemp1,
        &mut ftemp2,
        &r_params_int,
        &i_params_int,
    );

    // Evaluate at the midpoint
    i_params_int[2] = globals.em_steps as i32 + 2; // 1 = i_tm, 2 = k_x, and 3 = n_step

    i_params_int[3] = 0; // 0 -> s not equal to s_near or s_far; 10 -> s equal

    r_params_int[6] = del_x;

    func(
        x_curr + 0.5 * del_x,
        n_y,
        &mut f_out1,
        &mut f_out2,
        &r_params_int,
        &i_params_int,
    );

    // Combine everything into the evaluation of y_next
    for i in 0..n_y {
        y_curr[i] = y_curr[i] + f_out1[i] * del_x
            - f_out2[i] * del_x * (y_curr[i] + 0.5 * del_x * (ftemp1[i] - y_curr[i] * ftemp2[i]));
    }

    y_curr
}
