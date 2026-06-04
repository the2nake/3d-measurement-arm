#include <Eigen/Core>
#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numbers>
#include <print>
#include <random>
#include <stdexcept>

#include "GSA.hpp"
#include "GSA.ipp"

using namespace std;

class Eval {
 public:
  using eval_t = Eigen::Vector3d;
  using sense_t = Eigen::Vector3d;
  using param_t = Eigen::VectorXd;

  Eval(Eigen::Matrix3Xd targets,  // each column is a 3d target
       Eigen::Matrix3Xd senses,   // each column has 3 voltages
       std::function<eval_t(param_t, sense_t)> fk)
      : targets(targets), senses(senses), fk(fk) {
    if (senses.cols() != targets.cols()) {
      throw std::invalid_argument(
          std::format("dimension mismatch: senses.cols={}, targets.cols={}",
                      senses.cols(), targets.cols()));
    }
  }

  Eigen::Matrix3Xd effector_positions(const param_t& parameters) const {
    Eigen::Matrix3Xd positions(3, targets.cols());

    for (int i = 0; i < targets.cols(); ++i) {
      positions.col(i) = fk(parameters, senses.col(i));
    }
    return positions;
  }

  double operator()(const param_t& parameters) const {
    Eigen::Matrix3Xd residuals = targets - effector_positions(parameters);
    // use norm of mean absolute error
    return residuals.cwiseAbs().rowwise().sum().norm() / targets.cols();
  }

  const Eigen::Matrix3Xd targets;
  const Eigen::Matrix3Xd senses;
  const std::function<eval_t(param_t, sense_t)> fk;
};

constexpr double radians(double degrees) {
  return numbers::pi * degrees / 180.;
}

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
      {ct, -st * ca, st * sa,  a * ct},
      {st, ct * ca,  -ct * sa, a * st},
      {0., sa,       ca,       d     },
      {0., 0.,       0.,       1.    }
  };
}

Eigen::MatrixXd clean(Eigen::MatrixXd M) {
  M = (1e-6 < M.array().abs()).select(M, 0.);
  return M;
}

// param_vec_t = {origin_[x,y], d1, a2, a3, theta[1,2,3]_m, theta[1,2,3]_b}
//   where thetai = thetai_m * vi + thetai_b, vi is potentiometer value
// sense_vec_t = {v1, v2, v3}
using param_vec_t = Eigen::Vector<double, 11>;
using sense_vec_t = Eigen::Vector<double, 3>;

Eigen::Vector3d fk_gen(const param_vec_t& p, const sense_vec_t& v) {
  using namespace numbers;
  Eigen::Matrix4d T{
      {1., 0., 0., p(0)},
      {0., 1., 0., p(1)},
      {0., 0., 1., 0.  },
      {0., 0., 0., 1.  }
  };
  double j1 = p(5) * v(0) + p(8);
  double j2 = p(6) * v(1) + p(9);
  double j3 = p(7) * v(2) + p(10);
  // add pi / 2. to align the coordinate system
  T *= dh_matrix(0., p(2), pi / 2., pi / 2. + j1) *
       dh_matrix(p(3), 0., 0., j2) * dh_matrix(p(4), 0., 0., j3);
  return T(Eigen::seq(0, 2), 3);
};

int main() {
  using namespace numbers;

  // * note: test and calib datasets were collected with base local frame
  // rotated around +50 deg from test coordinate frame, with origin (40., -137.)
  // x, y parallel to ground, z is vertical
  const Eigen::Matrix<double, 3, 5> calib_points{
      {0., 100., 100., 0.,   50.},
      {0., 0.,   100., 100., 50.},
      {0., 0.,   0.,   0.,   0. }
  };
  const Eigen::Matrix<double, 3, 5> calib_volts{
      {5047., 6162., 5830., 5137., 5544.},
      {5688., 5686., 5792., 5775., 5710.},
      {4450., 4374., 3847., 3905., 4184.}
  };

  const Eigen::Matrix<double, 3, 5> test_points{
      {21., 35.5, 60.5, 77., 50. + 60. / sqrt(2)},
      {10., 43.5, 83.5, 40., 50. + 60. / sqrt(2)},
      {0.,  0.,   0.,   0.,  0.                 }
  };
  const Eigen::Matrix<double, 3, 5> test_volts{
      {5303., 5425., 5588., 5800., 5798.},
      {5687., 5702., 5751., 5703., 5775.},
      {4408., 4225., 3986., 4216., 3904.}
  };
  Eval eval_calib(calib_points, calib_volts, fk_gen);
  Eval eval_test(test_points, test_volts, fk_gen);

  // j1(v): 3704 -> 0,     6497 -> -pi/2 (fairly exact)
  // j2(v): 4777 -> pi/2,  6123 -> pi/4
  // j3(v): 2444 -> -pi/2, 3760 -> -3pi/4
  const double m1 = (-pi / 2.) / (6497 - 3704);
  const double m2 = (-pi / 4.) / (6123 - 4777);
  const double m3 = (-pi / 4.) / (3760 - 2444);
  const double b1 = 0 - m1 * 3704;
  const double b2 = pi / 2. - m2 * 4777;
  const double b3 = -pi / 2. - m3 * 2444;
  // * note: prior accurate up to variations in x, y, and b1
  param_vec_t prior = {
      30., -130., 62., 314., 333., m1, m2, m3, b1 + radians(50.), b2, b3};

  // print statistics for the initial estimate
  cout << endl;
  cout << "eval_calib(prior)=" << eval_calib(prior) << endl;
  cout << "residuals:\n"
       << calib_points - eval_calib.effector_positions(prior) << endl
       << endl;
  cout << "eval_test(prior)=" << eval_test(prior) << endl;
  cout << "residuals:\n"
       << test_points - eval_test.effector_positions(prior) << endl
       << endl;
  cout << endl;

  const double err_l = 5.;
  const double err_m = 1e-4;
  const double err_b = 0.1;
  param_vec_t dist = {15.,   15.,   err_l, err_l, err_l, err_m,
                      err_m, err_m, err_b, err_b, err_b};

  std::random_device r;
  // uncomment to seed using the random device
  // std::seed_seq ss{r(), r(), r(), r(), r(), r(), r(), r(), r()};
  // std::mt19937 rand(ss);
  std::mt19937 rand(1984);

  using Grav = GSA<param_vec_t, decltype(eval_calib)>;
  std::vector<param_vec_t> guesses =
      Grav::generate_guesses(prior, dist, rand, 200);
  Grav grav(guesses, eval_calib, rand);

  std::vector<param_vec_t> trace;
  string filename = "out/trace.txt";
  std::ofstream out_file(filename);
  if (!out_file.is_open()) {
    println("{} could not be opened", filename);
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();

  int best_iter = 0;
  double best_fitness = 1000.0;
  param_vec_t fit;
  while (grav.step()) {
    for (auto& pos : grav.positions()) { trace.emplace_back(pos); }

    // println();
    double curr_best = eval_calib(grav.best());
    if (curr_best < best_fitness) {
      best_fitness = curr_best;
      fit = grav.best();
      best_iter = grav.iter();
      println("[{}] new best, fitness = {:.3}", best_iter, best_fitness);
    }
    trace.emplace_back(grav.best());
  }

  auto end = std::chrono::high_resolution_clock::now();

  out_file << guesses.size() + 1 << endl;
  for (auto& pos : trace) { out_file << pos.transpose() << endl; }

  println();
  println("gsa: fitness = {:.3}, iter = {} / {}", best_fitness, best_iter,
          grav.m_max_iters);
  cout << "  score(test)=" << eval_test(fit) << endl;
  cout << "  params:" << endl << fit << endl;
  cout << "  residuals:" << endl
       << test_points - eval_test.effector_positions(fit) << endl;

  cout << endl
       << "optimisation took "
       << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
       << endl;
}
