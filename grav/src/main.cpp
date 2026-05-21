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

// todo: test calibeval, ideally on real case
class CalibEval {
 public:
  using eval_t = Eigen::Vector3d;
  using sense_t = Eigen::Vector3d;
  using param_t = Eigen::VectorXd;

  CalibEval(Eigen::Matrix3Xd targets,  // each column is a 3d target
            Eigen::Matrix3Xd senses,   // each column has 3 voltages
            std::function<eval_t(param_t, sense_t)> fk)
      : targets(targets), senses(senses), fk(fk) {
    if (senses.cols() != targets.cols()) {
      throw std::invalid_argument(
          std::format("dimension mismatch: senses.cols={}, targets.cols={}",
                      senses.cols(), targets.cols()));
    }
  }

  double operator()(const Eigen::MatrixXd& parameters) const {
    Eigen::Matrix3Xd effector_positions(3, targets.cols());

    for (int i = 0; i < targets.cols(); ++i) {
      effector_positions.col(i) = fk(parameters.col(i), senses.col(i));
    }

    Eigen::VectorXd residuals(targets.cols());
    for (int i = 0; i < targets.cols(); ++i) {
      residuals(i) = (targets.col(i) - effector_positions.col(i)).norm();
    }
    return residuals.norm();
  }

  const Eigen::Matrix3Xd targets;
  const Eigen::Matrix3Xd senses;
  const std::function<eval_t(param_t, sense_t)> fk;
};

/**
 * @brief Returns the homogeneous transform for Denavit-Hartenberg parameters
 *
 * This function uses a convention where the z-axis screw transform is applied
 * first
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

  // param_vec_t = {origin_[x,y], d1, a2, a3, theta[1,2,3]_m, theta[1,2,3]_b}
  // thetai = thetai_m * vi + thetai_b, vi is potentiometer value
  // sense_vec_t = {v1, v2, v3}
  using param_vec_t = Eigen::Vector<double, 11>;
  using sense_vec_t = Eigen::Vector<double, 3>;
  auto fk_gen = [](const param_vec_t& p,
                   const sense_vec_t& v) -> Eigen::Vector3d {
    Eigen::Matrix4d T{
        {1., 0., 0., p(0)},
        {0., 1., 0., p(1)},
        {0., 0., 1.,   0.},
        {0., 0., 0.,   1.}
    };
    double j1 = p(5) * v(0) + p(8);
    double j2 = p(6) * v(1) + p(9);
    double j3 = p(7) * v(2) + p(10);
    // add pi / 2. to align the coordinate system
    T *= dh_matrix(0., p(2), pi / 2., pi / 2. + j1) *
         dh_matrix(p(3), 0., 0., j2) * dh_matrix(p(4), 0., 0., j3);
    return T(Eigen::seq(0, 2), 3);
  };

  // x, y parallel to ground, z is vertical
  Eigen::Matrix<double, 3, 4> calib_points{
      {0., 100., 100.,   0.},
      {0.,   0., 100., 100.},
      {0.,   0.,   0.,   0.}
  };
  // todo: grab voltages from arm, maybe will need more inputs
  Eigen::Matrix<double, 3, 4> calib_volts{
      {0., 100., 100.,   0.},
      {0.,   0., 100., 100.},
      {0.,   0.,   0.,   0.}
  };
  CalibEval fk_eval(calib_points, calib_volts, fk_gen);

  // todo: get voltages and corresponding joint angles
  // todo: get calib_volts

  // j1(v): 3704 -> 0,     6497 -> -pi/2 (fairly exact)
  // j2(v): 4777 -> pi/2,  6123 -> pi/4
  // j3(v): 2444 -> -pi/2, 3760 -> -3pi/4
  double m1 = (-pi / 2.) / (6497 - 3704);
  double m2 = (-pi / 4.) / (6123 - 4777);
  double m3 = (-pi / 4.) / (3760 - 2444);
  double b1 = 0 - m1 * 3704;
  double b2 = pi / 2. - m2 * 4777;
  double b3 = -pi / 2. - m3 * 2444;
  param_vec_t mu_prior = {0., 0., 62., 314., 333., m1, m2, m3, b1, b2, b3};
  cout << endl;
  cout << "predicted location:" << endl
       << fk_gen(mu_prior, {5090., 5695., 4183.}) << endl;
  return 0;

  // todo: confirm that fk_eval(mu_prior) is low
  std::vector<param_vec_t> guesses = {mu_prior};

  GSA grav(guesses, fk_eval);
}

// ! fixme: sling issue where points diverge away from the center of mass

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
//   using namespace std::numbers;
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
