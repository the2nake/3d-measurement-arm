
#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <stdexcept>

#include "GSA.hpp"

using namespace std;

template <typename Vec, Evaluator<Vec> Eval>
  requires std::convertible_to<Vec, Eigen::VectorXd>
GSA<Vec, Eval>::GSA(const vector<Vec>& guesses, Eval err, std::mt19937 &rand)
    : m_eval(std::move(err)), m_x(guesses), m_v(guesses.size()), m_rand(rand) {
  if (guesses.size() < kb) {
    throw std::invalid_argument(
        std::format("not enough guesses, need at least {} (gsa)", kb));
  }
}

template <typename T, typename T2>
ostream& operator<<(ostream& os, const std::pair<T, T2>& p) {
  os << p.first << " " << p.second << " ";
  return os;
}

template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v) {
  for (const auto& x : v) { os << x << " "; }
  return os;
}

// ostream& operator<<(ostream& os, const vector<bool>& v) {
//   for (const auto& b : v) { os << (b ? "1" : "0") << " "; }
//   return os;
// }

template <typename Vec, Evaluator<Vec> Eval>
  requires std::convertible_to<Vec, Eigen::VectorXd>
bool GSA<Vec, Eval>::step() {
  // calculate error
  std::vector<double> fitnesses(m_x.size());
  double best_fitness = -1, worst_fitness = -1;
  for (int i = 0; i < m_x.size(); ++i) {
    double fitness = m_eval(m_x[i]);
    if (best_fitness == -1 || fitness < best_fitness) {
      best_fitness = fitness;
      m_best = m_x[i];
    }
    if (worst_fitness == -1 || fitness > worst_fitness) {
      worst_fitness = fitness;
    }
    fitnesses[i] = fitness;
  }
  assert(best_fitness != -1);
  assert(worst_fitness != -1);

  // compute non-normalised mass values
  double total = 0.0;
  std::vector<double> masses(fitnesses.size());
  for (int i = 0; i < fitnesses.size(); ++i) {
    masses[i] = (fitnesses[i] - worst_fitness) / (best_fitness - worst_fitness);
    total += masses[i];
  }

  // normalise mass values
  for (int i = 0; i < masses.size(); ++i) { masses[i] /= total; }

  // attract towards kb best scores
  std::vector<std::pair<double, int>> indexed_fitnesses(fitnesses.size());
  for (int i = 0; i < fitnesses.size(); ++i) {
    indexed_fitnesses[i] = {fitnesses[i], i};
  }
  std::nth_element(indexed_fitnesses.begin(),
                   indexed_fitnesses.begin() + kb - 1, indexed_fitnesses.end());

  // get gravitational constant
  const double G = G_i * std::exp(-beta * m_iter / (double)m_max_iters);
  std::vector<Vec> accels(m_x.size());
  std::uniform_real_distribution<> r(0, 1);
  for (int i = 0; i < m_x.size(); ++i) {
    for (int k = 0; k < kb; ++k) {
      int j = indexed_fitnesses[k].second;
      Vec d = m_x[j] - m_x[i];
      accels[i] +=
          r(m_rand) * G * masses[j] * (d) / (pow(d.norm(), rp) + epsilon);
    }
  }

  // TODO: use verlet integration?
  for (int i = 0; i < m_v.size(); ++i) { m_v[i] += accels[i]; m_x[i] += m_v[i]; }

  ++m_iter;
  return m_iter <= m_max_iters;
}
