/// Takes a cooled distribution of electrons and rebins them according to the
/// uncooled Lorentz factor array.
///
/// # Arguments
///
/// * `num_dp` - number of bins in the arrays
/// * `elec_gam_cooled` - Lorentz factor array for cooled distribution
/// * `elec_dist` - number of electrons (*not* dN/dp or dN/dgam) in each bin
/// * `elec_gam_array` - original Lorentz factor array
///
/// # Returns
///
/// * `elec_dist_rebin` - rebinned electron distribution
pub fn dist_rebin(
    num_dp: usize,
    elec_gam_cooled: &[f64],
    elec_dist: &[f64],
    elec_gam_array: &[f64],
) -> Vec<f64> {
    // Initialize variables
    let mut j_orig = 1usize;

    let mut elec_dist_rebin = vec![1.0e-99; num_dp];

    // The first bin is a gimme: since array entry 0 does not cool, the first
    // cooled bin is guaranteed to rebin entirely within the first bin of the
    // original Lorentz factor array
    if elec_dist[0] > 1.0e-55 {
        elec_dist_rebin[0] += elec_dist[0];
    }

    // Loop over bins in the cooled array and place them into the rebinned dist
    for i_cool in 1..num_dp {
        // Handle things differently if the current cooled bin falls entirely
        // within the uncooled bin under consideration.
        //
        // YES:
        //   Simply add the electrons to the rebinned distribution
        //
        // NO:
        //   Need to loop over bins in original Lorentz factor array until we
        //   have covered entirety of cooled array
        if elec_gam_cooled[i_cool] <= elec_gam_array[j_orig] {
            if elec_dist[i_cool] > 1.0e-55 {
                let elecs_to_add = elec_dist[i_cool];

                elec_dist_rebin[j_orig - 1] += elecs_to_add;
            }
        } else {
            while elec_gam_array[j_orig] <= elec_gam_cooled[i_cool] {
                let rebin_bottom = elec_gam_cooled[i_cool - 1].max(elec_gam_array[j_orig - 1]);

                let rebin_top = elec_gam_cooled[i_cool].min(elec_gam_array[j_orig]);

                // Compute number of electrons to place in new bin based on
                // fraction of elec_dist covered
                //
                //        |^^^^^^^^^^^^^^^^^^^|       elec_gam_cooled
                //  |_________|_________|_________|   elec_gam_array
                //         +++ +++++++++              <-- these parts
                if elec_dist[i_cool] > 1.0e-55 {
                    let elecs_to_add = elec_dist[i_cool] * (rebin_top - rebin_bottom)
                        / (elec_gam_cooled[i_cool] - elec_gam_cooled[i_cool - 1]);

                    elec_dist_rebin[j_orig - 1] += elecs_to_add;
                }

                j_orig += 1;
            }

            // Previous while-loop misses top part of elec_dist:
            //
            //        |^^^^^^^^^^^^^^^^^^^|       elec_gam_cooled
            //  |_________|_________|_________|   elec_gam_array
            //                       +++++ <-- this part
            //
            // To correct for that, add the missing electrons manually
            if elec_dist[i_cool] > 1.0e-55 {
                let elecs_to_add = elec_dist[i_cool]
                    * (elec_gam_cooled[i_cool] - elec_gam_array[j_orig - 1])
                    / (elec_gam_cooled[i_cool] - elec_gam_cooled[i_cool - 1]);

                elec_dist_rebin[j_orig - 1] += elecs_to_add;
            }
        }
    }

    elec_dist_rebin
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preserves_first_bin_when_uncooked() {
        let num_dp = 4;

        let elec_gam_cooled = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let elec_gam_array = vec![1.0, 2.0, 3.0, 4.0, 5.0];

        let elec_dist = vec![10.0, 0.0, 0.0, 0.0];

        let rebinned = dist_rebin(num_dp, &elec_gam_cooled, &elec_dist, &elec_gam_array);

        assert!((rebinned[0] - 10.0).abs() < 1.0e-12);
    }

    #[test]
    fn preserves_total_number_of_electrons() {
        let num_dp = 4;

        let elec_gam_cooled = vec![1.0, 2.2, 3.4, 4.6, 5.8];
        let elec_gam_array = vec![1.0, 2.0, 3.0, 4.0, 5.0];

        let elec_dist = vec![1.0, 2.0, 3.0, 4.0];

        let rebinned = dist_rebin(num_dp, &elec_gam_cooled, &elec_dist, &elec_gam_array);

        let original_sum: f64 = elec_dist.iter().sum();
        let rebinned_sum: f64 = rebinned.iter().sum();

        assert!(
            (original_sum - rebinned_sum).abs() < 1.0e-10,
            "original_sum = {original_sum}, rebinned_sum = {rebinned_sum}"
        );
    }

    #[test]
    fn rebinned_distribution_is_finite() {
        let num_dp = 4;

        let elec_gam_cooled = vec![1.0, 2.2, 3.4, 4.6, 5.8];
        let elec_gam_array = vec![1.0, 2.0, 3.0, 4.0, 5.0];

        let elec_dist = vec![1.0, 2.0, 3.0, 4.0];

        let rebinned = dist_rebin(num_dp, &elec_gam_cooled, &elec_dist, &elec_gam_array);

        assert!(rebinned.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn handles_empty_bins_correctly() {
        let num_dp = 4;

        let elec_gam_cooled = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let elec_gam_array = vec![1.0, 2.0, 3.0, 4.0, 5.0];

        let elec_dist = vec![0.0, 0.0, 5.0, 0.0];

        let rebinned = dist_rebin(num_dp, &elec_gam_cooled, &elec_dist, &elec_gam_array);

        let significant_bins = rebinned.iter().filter(|value| **value > 1.0e-55).count();

        assert_eq!(significant_bins, 1);
    }
}
