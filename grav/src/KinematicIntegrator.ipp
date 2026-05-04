#include <algorithm>

#include "KinematicIntegrator.hpp"

using namespace std;

template <NormedVector Vec, Evaluator<Vec> Eval>
GSA<Vec, Eval>::GSA(vector<Vec> points, Eval err, double dt)
    : dt(dt), err(err), m_prior_positions(points) {
  remove_duplicates(m_prior_positions);
  m_positions = vector<Vec>(m_prior_positions.size());
  m_masses = vector<double>(m_prior_positions.size());
  m_accels = vector<Vec>(m_prior_positions.size());
  compute_forces(m_prior_positions);
  for (int i = 0; i < m_prior_positions.size(); ++i) {
    m_positions[i] = m_prior_positions[i] + 0.5 * m_accels[i] * dt * dt;
  }
}

template <NormedVector Vec, Evaluator<Vec> Eval>
void GSA<Vec, Eval>::remove_duplicates(vector<Vec>& positions) {
  vector<Vec> new_positions;
  new_positions.reserve(positions.size());
  vector<bool> enabled(positions.size(), true);

  for (int i = 0; i < positions.size(); ++i) {
    if (!enabled[i]) continue;
    for (int j = i + 1; j < positions.size(); ++j) {
      if ((positions[i] - positions[j]).norm() < min_dist) {
        enabled[i] = false;
      }
    }
  }

  for (int i = 0; i < positions.size(); ++i) {
    if (!enabled[i]) continue;
    new_positions.emplace_back(positions[i]);
  }
  positions = new_positions;
}

template <NormedVector Vec, Evaluator<Vec> Eval>
void GSA<Vec, Eval>::compute_forces(const vector<Vec>& positions) {
  for (auto& force : m_accels) { force = 0.0 * force; }
  for (int i = 0; i < positions.size(); ++i) {
    m_masses[i] = 1.0 / err(positions[i]);
  }

  // TODO: architecture change, compute_forces should just operate on
  // m_positions and m_prior_positions directly should integrate this bit into
  // force computation? would this lead to invalidation, or figure out how to
  // compensate for impact on previously computed pairs. should be fine to
  // merge and keep going, though.

  for (int i = 0; i < positions.size(); ++i) {
    // TODO: mark a "skip" flag for deleted objects
    for (int j = i + 1; j < positions.size(); ++j) {
      double distance = (positions[i] - positions[j]).norm();
      // if (distance < min_dist) {
      //   // TODO: architecture change affects this bit of code
      //   positions[i] =
      //       (m_masses[i] * positions[i] + m_masses[j] * positions[j]) /
      //       (m_masses[i] + m_masses[j]);
      //   // conservation of momentum
      //   m_prior_positions[i] = (m_masses[i] * m_prior_positions[i] +
      //                           m_masses[j] * m_prior_positions[j]) /
      //                          (m_masses[i] + m_masses[j]);
      //   // TODO: mark object at index `j` for deletiondd
      // }
      // TODO! remove the need for the hacky ==> if (distance < 0.1) continue;
      double magnitude = m_masses[i] * m_masses[j] / distance;
      Vec force_i_to_j = magnitude * (positions[j] - positions[i]) / distance;
      m_accels[i] += force_i_to_j;
      m_accels[j] -= force_i_to_j;
    }
  }

  for (int i = 0; i < positions.size(); ++i) {
    m_accels[i] = m_accels[i] / m_masses[i];
  }
}

template <NormedVector Vec, Evaluator<Vec> Eval>
void GSA<Vec, Eval>::step() {
  compute_forces(m_positions);

  // using Verlet integration
  const auto starting_positions = m_positions;

  for (int i = 0; i < m_positions.size(); ++i) {
    // TODO: add support for variable simulation timestep based on maximuum
    // acceleration value
    m_positions[i] =
        2 * m_positions[i] - m_prior_positions[i] + m_accels[i] * dt * dt;
  }
  m_prior_positions = starting_positions;
}

template <NormedVector Vec, Evaluator<Vec> Eval>
Vec GSA<Vec, Eval>::best() const {
  return *min_element(m_positions.begin(), m_positions.end(),
                      [this](auto& a, auto& b) { return err(a) < err(b); });
}
