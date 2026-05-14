#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <format>
#include <fstream>
#include <functional>
#include <numbers>
#include <print>
#include <random>
#include <stdexcept>
#include <type_traits>

#include "GSA.hpp"
#include "GSA.ipp"

using namespace std;

class CalibEval {
 public:
  using eval_t = Eigen::Vector3d;
  using param_t = Eigen::VectorXd;

  CalibEval(Eigen::Matrix3Xd targets,  // each column is a 3d target
            std::function<eval_t(param_t)> fk)
      : targets(targets), fk(fk) {}

  double operator()(const Eigen::MatrixXd& parameters) const {
    if (parameters.cols() != targets.cols())
      throw std::invalid_argument(
          std::format("incorrect input dimension: got {}, expected {}",
                      parameters.cols(), targets.cols()));
    Eigen::Matrix3Xd effector_positions(3, targets.cols());

    for (int i = 0; i < targets.cols(); ++i) {
      effector_positions.col(i) = fk(parameters.col(i));
    }

    Eigen::VectorXd residuals(targets.cols());
    for (int i = 0; i < targets.cols(); ++i) {
      residuals(i) = (targets.col(i) - effector_positions.col(i)).norm();
    }
    return residuals.norm();
  }

  const Eigen::Matrix3Xd targets;
  const std::function<eval_t(param_t)> fk;
};

/**
 * @brief Returns the homogeneous transform for Denavit-Hartenberg parameters
 *
 * @param a link length
 * @param d link offset
 * @param alpha link twist in radians
 * @param theta joint angle in radians
 * @return Eigen::Matrix4d containing the transform. Transforms should be
 * chained in increasing order from left to right (e.g. `A_1_2 = A_1 * A_2`)
 */
Eigen::Matrix4d dh_matrix(double a, double d, double alpha, double theta) {
  double ct = std::cos(theta);
  double st = std::sin(theta);

  double ca = std::cos(alpha);
  double sa = std::sin(alpha);

  return Eigen::Matrix4d{
      {ct, -st * ca,  st * sa, a * ct},
      {st,  ct * ca, -ct * sa, a * st},
      {0.,       sa,       ca,      d},
      {0.,       0.,       0.,     1.}
  };
}
Eigen::MatrixXd clean(Eigen::MatrixXd M) {
  M = (1e-6 < M.array().abs()).select(M, 0.);
  return M;
}

int main() {
  using namespace numbers;

  double j1 = 0.;
  double j2 = 45. * pi / 180.;
  double j3 = -45. * pi / 180.;
  // 90 deg offset is important
  Eigen::MatrixXd A1 = dh_matrix(0., 10., pi / 2., pi / 2 + j1);
  Eigen::MatrixXd A2 = dh_matrix(10., 0., 0., j2);
  Eigen::MatrixXd A3 = dh_matrix(10., 0., 0., j3);
  cout << clean(A1 * A2 * A3) << endl;
  return 0;

  using param_vec_t = Eigen::Vector<double, 6>;
  // todo: implmenet fk_gen. should compute chain of transformations using
  // dh_matrix and return end-effector coordinates (T.col(3).)
  auto fk_gen = [](const param_vec_t& v) {
    return Eigen::Vector3d{0.0, 0.0, 0.0};
  };
  // todo: make composite evaluator, taking in matrices and evaluating each
  // column, operator() overloaded as norm of vector of errors
  // x, y parallel to ground, z is vertical
  Eigen::Matrix<double, 3, 4> calibration_points{
      {0., 100., 100.,   0.},
      {0.,   0., 100., 100.},
      {0.,   0.,   0.,   0.}
  };
  CalibEval fk_eval(calibration_points, fk_gen);

  param_vec_t mu_prior = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<param_vec_t> guesses = {mu_prior};

  GSA grav(guesses, fk_eval);
}

// sling issue where points diverge away from the center of mass
// todo: email the people who originally did the GSA method and ask how they
// handled the slingshot issue

// template <typename TargetT, typename SpaceT>
// class LstsqEval {
//  public:
//   LstsqEval(TargetT target, std::function<TargetT(SpaceT)> gen)
//       : target(target), gen(gen) {}
//   double operator()(SpaceT v) const {
//     TargetT output = gen(v);
//     double norm = (target - output).norm();
//     return norm;
//   };
//   const TargetT target;
//   const std::function<TargetT(SpaceT)> gen;
// };

// std::vector<Eigen::Vector2d> populate_initial_guesses(int samples = 9,
//                                                       double var = 0.01) {
//   Eigen::Vector2d target(3, 3);
//   const double E_mod = 3 * sqrt(2), V_mod = var * E_mod;
//   const double E_arg = pi / 4, V_arg = var * E_arg;
//   std::random_device r;
//   std::seed_seq ss{r(), r(), r(), r(), r(), r(), r(), r(), r()};
//   std::mt19937 mt(ss);
//   std::uniform_real_distribution<> dist_mod(E_mod - V_mod, E_mod + V_mod);
//   std::uniform_real_distribution<> dist_arg(E_arg - V_arg, E_arg + V_arg);
//   std::vector<Eigen::Vector2d> guesses;
//   for (int i = 0; i < samples; ++i) {
//     guesses.emplace_back(dist_mod(mt), dist_arg(mt));
//   }
//   // guesses.emplace_back(E_mod, E_arg);
//   // guesses.emplace_back(E_mod, E_arg);
//   return guesses;
// }

// int main() {
//   Eigen::Vector2d target(3, 3);
//   std::vector<Eigen::Vector2d> guesses = populate_initial_guesses(150, 0.01);

//   auto generator = [](const Eigen::Vector2d& v2) -> Eigen::Vector2d {
//     return {v2[0] * std::cos(v2[1]), v2[0] * std::sin(v2[1])};
//   };
//   LstsqEval<Eigen::Vector2d, Eigen::Vector2d> eval(target, generator);

//   GSA grav(guesses, eval);

//   std::vector<Eigen::Vector2d> trace;
//   string filename = "out/trace.txt";
//   std::ofstream out_file(filename);
//   if (!out_file.is_open()) {
//     println("{} could not be opened", filename);
//     return 1;
//   }

//   int best_iter = 0;
//   double best_score = 1000.0;
//   Eigen::Vector2d fit;
//   while (grav.step()) {
//     for (auto& pos : grav.positions()) { trace.emplace_back(pos); }

//     // println();
//     double new_best = eval(grav.best());
//     if (new_best < best_score) {
//       best_score = new_best;
//       fit = grav.best();
//       best_iter = grav.iter();
//       println("[{}] new best, err = {:.3}", best_iter, best_score);
//     }
//     trace.emplace_back(grav.best());
//   }

//   out_file << guesses.size() + 1 << endl;
//   for (auto& pos : trace) { out_file << pos[0] << " " << pos[1] << endl; }

//   println();
//   println("    gsa: params: {{{: 3.5f}, {: 3.5f}}}, err = {:.3}, iter = {}",
//           fit[0], fit[1], best_score, best_iter);
//   println("optimal: params: {{{: 3.5f}, {: 3.5f}}}", sqrt(18.), atan2(3, 3));

//   system("python3 scripts/plotter.py");
// }
