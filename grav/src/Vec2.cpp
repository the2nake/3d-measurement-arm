#include "Vec2.hpp"

bool operator==(const Vec2& a, const Vec2& b) {
  return a.x1 == b.x1 && a.x2 == b.x2;
}
Vec2 operator*(double x, const Vec2& a) { return {x * a.x1, x * a.x2}; }
Vec2 operator*(const Vec2& a, double x) { return operator*(x, a); }
Vec2 operator/(const Vec2& a, double x) { return (1 / x) * a; }
Vec2 operator+(const Vec2& a, const Vec2& b) {
  return {a.x1 + b.x1, a.x2 + b.x2};
}
Vec2 operator-(const Vec2& a, const Vec2& b) { return a + (-1.0 * b); }
