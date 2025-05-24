#include <RcppEigen.h>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List omp_encode_rcpp(
    const Eigen::Map<Eigen::VectorXd> signal_y,
    const Eigen::MappedSparseMatrix<double> dict_D,
    double residual_norm_sq_tol,
    int max_active_atoms_L)
{
    Eigen::VectorXd residual = signal_y;
    std::vector<int> active_indices;
    Eigen::VectorXd coeffs;

    int iter = 0;
    while (residual.squaredNorm() > residual_norm_sq_tol &&
           iter < max_active_atoms_L) {
        // Correlation with dictionary atoms
        Eigen::VectorXd corr = dict_D.transpose() * residual;
        Eigen::Index best_j;
        corr.cwiseAbs().maxCoeff(&best_j);

        if (std::find(active_indices.begin(), active_indices.end(), best_j) == active_indices.end()) {
            active_indices.push_back(static_cast<int>(best_j));
        }

        Eigen::SparseMatrix<double> Dsub(dict_D.rows(), active_indices.size());
        std::vector<Eigen::Triplet<double>> trips;
        for (size_t col = 0; col < active_indices.size(); ++col) {
            int j = active_indices[col];
            for (Eigen::MappedSparseMatrix<double>::InnerIterator it(dict_D, j); it; ++it) {
                trips.emplace_back(it.row(), static_cast<int>(col), it.value());
            }
        }
        if (!trips.empty()) {
            Dsub.setFromTriplets(trips.begin(), trips.end());
        }

        Eigen::SparseMatrix<double> gram = Dsub.transpose() * Dsub;
        Eigen::VectorXd rhs = Dsub.transpose() * signal_y;

        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(gram);
        if (solver.info() != Eigen::Success) {
            stop("SimplicialLLT decomposition failed");
        }
        coeffs = solver.solve(rhs);
        residual = signal_y - Dsub * coeffs;
        ++iter;
    }

    bool max_iter_reached = (iter == max_active_atoms_L) &&
                            (residual.squaredNorm() > residual_norm_sq_tol);
    std::vector<double> coeff_vec(coeffs.data(), coeffs.data() + coeffs.size());

    return Rcpp::List::create(
        Rcpp::Named("indices_0based") = active_indices,
        Rcpp::Named("coefficients") = coeff_vec,
        Rcpp::Named("max_iter_reached") = max_iter_reached);
}

