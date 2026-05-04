#pragma once
#include <concepts>
#include <vector>

template <typename T>
concept NormedVector = requires(T a, T b, double scalar) {
  0.0 * a == 0.0 * b;
  0.0 * a == T();
  std::equality_comparable<T>;
  { a + b } -> std::same_as<T>;
  { a - b } -> std::same_as<T>;
  { scalar * a } -> std::same_as<T>;
  { a.norm() } -> std::convertible_to<double>;
};

template <typename T, typename Vec>
concept Evaluator = requires(T a, Vec v) {
  { a(v) } -> std::convertible_to<double>;
};

template <NormedVector Vec, Evaluator<Vec> Eval>
class GSA {
 public:
  GSA(std::vector<Vec> points, Eval metric, double dt);

  void step();
  Vec best() const;

  const double G = 0.5;
  const double dt;
  const Eval err;

  const std::vector<Vec>& positions() const { return m_positions; }

 private:
  void remove_duplicates(std::vector<Vec>& positions);
  void compute_forces(const std::vector<Vec>& positions);

  std::vector<Vec> m_positions;
  std::vector<Vec> m_prior_positions;
  std::vector<Vec> m_accels;
  std::vector<double> m_masses;
};
