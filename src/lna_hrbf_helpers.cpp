#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix poisson_disk_sample_component_rcpp(IntegerMatrix component_vox_coords_0based, double radius_vox_sq, int component_seed) {
    RNGScope scope; // ensure RNG gets properly set
    Function set_seed("set.seed");
    set_seed(component_seed);

    int nvox = component_vox_coords_0based.nrow();
    if (nvox == 0) {
        return IntegerMatrix(0, 3);
    }

    // Create permutation of indices using R's sample.int
    Function sample_int("sample.int");
    IntegerVector perm = sample_int(nvox, nvox, false);

    std::vector< std::array<int,3> > selected;
    selected.reserve(nvox); // worst case

    for (int idx = 0; idx < nvox; ++idx) {
        int row = perm[idx] - 1; // convert 1-based to 0-based index
        int i = component_vox_coords_0based(row, 0);
        int j = component_vox_coords_0based(row, 1);
        int k = component_vox_coords_0based(row, 2);

        bool accept = true;
        for (const auto& centre : selected) {
            double dx = static_cast<double>(i - centre[0]);
            double dy = static_cast<double>(j - centre[1]);
            double dz = static_cast<double>(k - centre[2]);
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < radius_vox_sq) {
                accept = false;
                break; // early exit
            }
        }
        if (accept) {
            selected.push_back({i,j,k});
        }
    }

    int nsel = selected.size();
    IntegerMatrix out(nsel, 3);
    for (int r = 0; r < nsel; ++r) {
        out(r,0) = selected[r][0];
        out(r,1) = selected[r][1];
        out(r,2) = selected[r][2];
    }
    return out;
}

