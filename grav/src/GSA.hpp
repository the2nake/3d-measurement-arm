#pragma once
#include <concepts>
#include <random>
#include <vector>

template <typename T>
concept NormedVector = requires(T a, T b, double scalar) {
  std::equality_comparable<T>;
  { a + b } -> std::same_as<T>;
  { a - b } -> std::same_as<T>;
  { scalar * a } -> std::same_as<T>;
  { a.norm() } -> std::convertible_to<double>;
  0.0 * a == 0.0 * b;
  0.0 * a == T();  // should default to zero vector
};

template <typename T, typename Vec>
concept Evaluator = requires(T a, Vec v) {
  a(v) >= 0.0;
  { a(v) } -> std::convertible_to<double>;
};

template <NormedVector Vec, Evaluator<Vec> Eval>
class GSA {
 public:
  GSA(std::vector<Vec> guesses, Eval metric);

  bool step();
  Vec best() const;

  int iter() const { return m_iter; }
  const std::vector<Vec>& positions() const { return m_x; }

  const Eval m_eval;

  const int m_max_iters = 1e4;
  const int rp = 1;         // exponent f euclidean distance
  const int kb = 2;         // number of best solutions to pick
  const double G_i = 1e-5;  // initial gravitation, how to choose this?
  const double beta = 5.0;  // increasing beta increases gravitation falloff
  const double epsilon = 0.000001;

 private:
  std::mt19937 m_rand;

  int m_iter = 0;

  std::vector<Vec> m_x;  // positions
  std::vector<Vec> m_v;  // velocities
};
