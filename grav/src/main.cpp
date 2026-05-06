#include <cmath>
#include <fstream>
#include <functional>
#include <numbers>
#include <print>
#include <random>

#include "GSA.hpp"
#include "GSA.ipp"
#include "Vec2.hpp"

using namespace std;

// sling issue where points diverge away from the center of mass

template <NormedVector TargetT, NormedVector SpaceT>
class ModArgEval {
 public:
  ModArgEval(TargetT target, std::function<TargetT(SpaceT)> gen)
      : target(target), gen(gen) {}
  double operator()(SpaceT v) const {
    TargetT output = gen(v);
    double norm = (target - output).norm();
    return norm * norm;
  };
  const TargetT target;
  const std::function<TargetT(SpaceT)> gen;
};

std::vector<Vec2> populate_initial_guesses(int samples = 9, double var = 0.01) {
  Vec2 target(3, 3);
  const double E_mod = 3 * sqrt(2), V_mod = var * E_mod;
  const double E_arg = numbers::pi / 4, V_arg = var * E_arg;
  std::random_device r;
  std::seed_seq ss{r(), r(), r(), r(), r(), r(), r(), r(), r()};
  std::mt19937 mt(ss);
  std::uniform_real_distribution<> dist_mod(E_mod - V_mod, E_mod + V_mod);
  std::uniform_real_distribution<> dist_arg(E_arg - V_arg, E_arg + V_arg);
  std::vector<Vec2> guesses;
  for (int i = 0; i < samples; ++i) {
    guesses.emplace_back(dist_mod(mt), dist_arg(mt));
  }
  // guesses.emplace_back(E_mod, E_arg);
  // guesses.emplace_back(E_mod, E_arg);
  return guesses;
}

int main() {
  Vec2 target(3, 3);

  std::vector<Vec2> guesses = populate_initial_guesses(150, 0.1);

  auto generator = [](const Vec2& v2) {
    return Vec2{v2.x1 * std::cos(v2.x2), v2.x1 * std::sin(v2.x2)};
  };
  ModArgEval<Vec2, Vec2> eval(target, generator);
  GSA<Vec2, ModArgEval<Vec2, Vec2>> grav(guesses, eval);

  std::vector<Vec2> trace;
  string filename = "out/trace.txt";
  std::ofstream out_file(filename);
  if (!out_file.is_open()) {
    println("{} could not be opened", filename);
    return 1;
  }

  int best_iter = 0;
  double best_score = 1000.0;
  Vec2 fit;
  while (grav.step()) {
    for (auto& pos : grav.positions()) { trace.emplace_back(pos); }

    // println();
    double new_best = eval(grav.best());
    if (new_best < best_score) {
      best_score = new_best;
      fit = grav.best();
      best_iter = grav.iter();
      println("[{}] new best, err^2 = {:.3}", best_iter, best_score);
    }
    trace.emplace_back(grav.best());
  }

  out_file << guesses.size() + 1 << endl;
  for (auto& pos : trace) { out_file << pos.x1 << " " << pos.x2 << endl; }

  println();
  println("    gsa: params: {{{: 3.5f}, {: 3.5f}}}, err^2 = {:.3}, iter = {}",
          fit.x1, fit.x2, best_score, best_iter);
  println("optimal: params: {{{: 3.5f}, {: 3.5f}}}", sqrt(18.), atan2(3, 3));

  system("python3 scripts/plotter.py");
}
