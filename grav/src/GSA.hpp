#pragma once
#include <Eigen/Dense>
#include <concepts>
#include <random>
#include <vector>

template <typename T, typename Vec>
concept Evaluator = requires(T a, Vec v) {
  a(v) >= 0.0;
  { a(v) } -> std::convertible_to<double>;
};

template <typename Vec, Evaluator<Vec> Eval>
  requires std::convertible_to<Vec, Eigen::VectorXd>
class GSA {
 public:
  template <typename T>
  static std::vector<T> generate_guesses(T mean, T range, std::mt19937& rand,
                                         int n = 150) {
    if (mean.rows() != range.rows()) {
      throw std::invalid_argument("mismatched input dimensions");
    }
    std::uniform_real_distribution<> uni(-0.5, 0.5);  // use this then scaled

    std::vector<T> guesses = {mean};

    for (int i = 0; i < n; ++i) {
      T guess(mean.rows());
      for (int j = 0; j < mean.rows(); ++j) {
        guess(j) = mean(j) + uni(rand) * range(j);
      }
      guesses.emplace_back(guess);
    }

    return guesses;
  }

  GSA(const std::vector<Vec>& guesses, Eval metric, std::mt19937& rand);

  bool step();

  int iter() const { return m_iter; }
  Vec best() const { return m_best; }
  const std::vector<Vec>& positions() const { return m_x; }

  const Eval m_eval;

  const int m_max_iters = 1e3;   // low impact
  const int rp = 1;              // exponent of euclidean distance, low impact
  const int kb = 2;              // number of best solutions to pick, low impact
  const double G_i = 1e-1;       // how to choose?, high impact
  const double beta = 20.0;      // promotes gravitation falloff, medium impact
  const double epsilon = 1e-60;  // lower is better, medium impact

 private:
  double pow(double x, int p) { return p == 0 ? 1 : x * pow(x, p - 1); }

  std::mt19937& m_rand;

  int m_iter = 0;
  Vec m_best;

  std::vector<Vec> m_x;  // positions
  std::vector<Vec> m_v;  // velocities
};
