#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::SparseMatrix<float> hrbf_atoms_rcpp(
    const Eigen::Map<Eigen::MatrixXf> mask_xyz_world,
    const Eigen::Map<Eigen::MatrixXf> centres_xyz_world,
    const Eigen::Map<Eigen::VectorXf> sigma_vec_mm,
    std::string kernel_type,
    double value_threshold = 1e-8)
{
    const int N = mask_xyz_world.rows();
    const int K = centres_xyz_world.rows();
    if (sigma_vec_mm.size() != K) {
        Rcpp::stop("sigma_vec_mm length must match centres rows");
    }
    const bool gaussian = (kernel_type == "gaussian");

    std::vector< Eigen::Triplet<float> > triplets;
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    std::vector< std::vector<Eigen::Triplet<float>> > thread_triplets(nthreads);
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::vector<Eigen::Triplet<float>>& local = thread_triplets[tid];
#pragma omp for schedule(static)
        for (int k = 0; k < K; ++k) {
            float sigma = sigma_vec_mm[k];
            float inv_two_sigma2 = 1.0f / (2.0f * sigma * sigma);
            float inv_sigma = 1.0f / sigma;
            for (int n = 0; n < N; ++n) {
                float dx = mask_xyz_world(n,0) - centres_xyz_world(k,0);
                float dy = mask_xyz_world(n,1) - centres_xyz_world(k,1);
                float dz = mask_xyz_world(n,2) - centres_xyz_world(k,2);
                float dist2 = dx*dx + dy*dy + dz*dz;
                float value;
                if (gaussian) {
                    value = std::exp(-dist2 * inv_two_sigma2);
                } else {
                    float dist = std::sqrt(dist2) * inv_sigma;
                    if (dist >= 1.0f) continue;
                    float base = 1.0f - dist;
                    value = std::pow(base,8.0f) * (32.0f*dist*dist*dist +
                                                  25.0f*dist*dist + 8.0f*dist + 1.0f);
                }
                if (std::fabs(value) > value_threshold) {
                    local.emplace_back(k, n, value);
                }
            }
        }
    }
    for (auto &vec : thread_triplets) {
        triplets.insert(triplets.end(), vec.begin(), vec.end());
    }
#else
    triplets.reserve(static_cast<size_t>(N) * 4);
    for (int k = 0; k < K; ++k) {
        float sigma = sigma_vec_mm[k];
        float inv_two_sigma2 = 1.0f / (2.0f * sigma * sigma);
        float inv_sigma = 1.0f / sigma;
        for (int n = 0; n < N; ++n) {
            float dx = mask_xyz_world(n,0) - centres_xyz_world(k,0);
            float dy = mask_xyz_world(n,1) - centres_xyz_world(k,1);
            float dz = mask_xyz_world(n,2) - centres_xyz_world(k,2);
            float dist2 = dx*dx + dy*dy + dz*dz;
            float value;
            if (gaussian) {
                value = std::exp(-dist2 * inv_two_sigma2);
            } else {
                float dist = std::sqrt(dist2) * inv_sigma;
                if (dist >= 1.0f) continue;
                float base = 1.0f - dist;
                value = std::pow(base,8.0f) * (32.0f*dist*dist*dist +
                                              25.0f*dist*dist + 8.0f*dist + 1.0f);
            }
            if (std::fabs(value) > value_threshold) {
                triplets.emplace_back(k, n, value);
            }
        }
    }
#endif
    Eigen::SparseMatrix<float> result(K, N);
    if (!triplets.empty()) {
        result.setFromTriplets(triplets.begin(), triplets.end());
    }
    return result;
}
