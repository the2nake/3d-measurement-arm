#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <functional>
#include <numbers>
#include <print>
#include <random>

#include "GSA.hpp"
#include "GSA.ipp"

using namespace std;

// sling issue where points diverge away from the center of mass

template <typename TargetT, typename SpaceT>
class LstsqEval {
 public:
  LstsqEval(TargetT target, std::function<TargetT(SpaceT)> gen)
      : target(target), gen(gen) {}
  double operator()(SpaceT v) const {
    TargetT output = gen(v);
    double norm = (target - output).norm();
    return norm;
  };
  const TargetT target;
  const std::function<TargetT(SpaceT)> gen;
};

std::vector<Eigen::Vector2d> populate_initial_guesses(int samples = 9,
                                                      double var = 0.01) {
  Eigen::Vector2d target(3, 3);
  const double E_mod = 3 * sqrt(2), V_mod = var * E_mod;
  const double E_arg = numbers::pi / 4, V_arg = var * E_arg;
  std::random_device r;
  std::seed_seq ss{r(), r(), r(), r(), r(), r(), r(), r(), r()};
  std::mt19937 mt(ss);
  std::uniform_real_distribution<> dist_mod(E_mod - V_mod, E_mod + V_mod);
  std::uniform_real_distribution<> dist_arg(E_arg - V_arg, E_arg + V_arg);
  std::vector<Eigen::Vector2d> guesses;
  for (int i = 0; i < samples; ++i) {
    guesses.emplace_back(dist_mod(mt), dist_arg(mt));
  }
  // guesses.emplace_back(E_mod, E_arg);
  // guesses.emplace_back(E_mod, E_arg);
  return guesses;
}

int main() {
  Eigen::Vector2d target(3, 3);
  std::vector<Eigen::Vector2d> guesses = populate_initial_guesses(150, 0.01);

  auto generator = [](const Eigen::Vector2d& v2) -> Eigen::Vector2d {
    return {v2[0] * std::cos(v2[1]), v2[0] * std::sin(v2[1])};
  };
  LstsqEval<Eigen::Vector2d, Eigen::Vector2d> eval(target, generator);

  GSA grav(guesses, eval);

  std::vector<Eigen::Vector2d> trace;
  string filename = "out/trace.txt";
  std::ofstream out_file(filename);
  if (!out_file.is_open()) {
    println("{} could not be opened", filename);
    return 1;
  }

  int best_iter = 0;
  double best_score = 1000.0;
  Eigen::Vector2d fit;
  while (grav.step()) {
    for (auto& pos : grav.positions()) { trace.emplace_back(pos); }

    // println();
    double new_best = eval(grav.best());
    if (new_best < best_score) {
      best_score = new_best;
      fit = grav.best();
      best_iter = grav.iter();
      println("[{}] new best, err = {:.3}", best_iter, best_score);
    }
    trace.emplace_back(grav.best());
  }

  out_file << guesses.size() + 1 << endl;
  for (auto& pos : trace) { out_file << pos[0] << " " << pos[1] << endl; }

  println();
  println("    gsa: params: {{{: 3.5f}, {: 3.5f}}}, err = {:.3}, iter = {}",
          fit[0], fit[1], best_score, best_iter);
  println("optimal: params: {{{: 3.5f}, {: 3.5f}}}", sqrt(18.), atan2(3, 3));

  system("python3 scripts/plotter.py");
}
