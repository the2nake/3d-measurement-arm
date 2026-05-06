#pragma once

struct Vec2 {
  Vec2() : x1(0.0), x2(0.0) {}
  Vec2(double mod, double arg) : x1(mod), x2(arg) {}

  double x1;
  double x2;

  double norm() { return x1 * x1 + x2 * x2; }

  void operator+=(const Vec2& b) {
    x1 += b.x1;
    x2 += b.x2;
  }

  void operator-=(const Vec2& b) {
    x1 -= b.x1;
    x2 -= b.x2;
  }
};

bool operator==(const Vec2& a, const Vec2& b);
Vec2 operator*(double x, const Vec2& a);
Vec2 operator*(const Vec2& a, double x);
Vec2 operator/(const Vec2& a, double x);
Vec2 operator+(const Vec2& a, const Vec2& b);
Vec2 operator-(const Vec2& a, const Vec2& b);
