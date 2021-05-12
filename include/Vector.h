#pragma once
#include <array>
#include <cmath>
#include <iostream>

namespace dem
{
template <class Type>
struct Vector2D
{
  template <class U>
  using vector_type = Vector2D<U>;
  using value_type = Type;
  value_type x, y;

  Vector2D() noexcept = default;
  constexpr Vector2D(const Vector2D&) noexcept = default;

  template <class X, class Y>
  constexpr Vector2D(X _X, Y _Y) noexcept
    : x(static_cast<value_type>(_X)), y(static_cast<value_type>(_Y))
  {
  }

  constexpr Vector2D(value_type _X, value_type _Y) noexcept : x(_X), y(_Y) {}

  template <class U>
  constexpr Vector2D(const Vector2D<U>& V) noexcept
    : x(static_cast<value_type>(V.x)), y(static_cast<value_type>(V.y))
  {
  }

  constexpr Vector2D& Set(value_type _X, value_type _Y) noexcept
  {
    x = _X;
    y = _Y;
    return *this;
  }

  constexpr Vector2D& Set(const Vector2D& V) noexcept { return *this = V; }

  constexpr Vector2D operator+() const noexcept { return *this; }

  constexpr Vector2D operator-() const noexcept { return {-x, -y}; }

  constexpr Vector2D operator+(const Vector2D& Other) const noexcept
  {
    return {x + Other.x, y + Other.y};
  }

  constexpr Vector2D operator-(const Vector2D& Other) const noexcept
  {
    return {x - Other.x, y - Other.y};
  }

  constexpr Vector2D operator*(const Vector2D& Other) const noexcept
  {
    return {x * Other.x, y * Other.y};
  }

  constexpr Vector2D operator/(const Vector2D& Other) const noexcept
  {
    return {x / Other.x, y / Other.y};
  }

  constexpr Vector2D operator+(const std::array<int, 2>& Other) const noexcept
  {
    return {x + Other[0], y + Other[1]};
  }

  constexpr Vector2D operator-(const std::array<int, 2>& Other) const noexcept
  {
    return {x - Other[0], y - Other[1]};
  }

  constexpr Vector2D operator*(const std::array<int, 2>& Other) const noexcept
  {
    return {x * Other[0], y * Other[1]};
  }

  constexpr Vector2D operator/(const std::array<int, 2>& Other) const noexcept
  {
    return {x / Other[0], y / Other[1]};
  }

  constexpr Vector2D operator+(value_type S) const noexcept { return {x + S, y + S}; }

  constexpr Vector2D operator-(value_type S) const noexcept { return {x - S, y - S}; }

  constexpr Vector2D operator*(value_type S) const noexcept { return {x * S, y * S}; }

  constexpr Vector2D operator/(value_type S) const noexcept
  {
    return *this * (static_cast<value_type>(1.0) / S);
  }

  constexpr Vector2D& operator+=(const Vector2D& Other) noexcept
  {
    x += Other.x;
    y += Other.y;
    return *this;
  }

  constexpr Vector2D& operator-=(const Vector2D& Other) noexcept
  {
    x -= Other.x;
    y -= Other.y;
    return *this;
  }

  constexpr Vector2D& operator*=(const Vector2D& Other) noexcept
  {
    x *= Other.x;
    y *= Other.y;
    return *this;
  }

  constexpr Vector2D& operator/=(const Vector2D& Other) noexcept
  {
    x /= Other.x;
    y /= Other.y;
    return *this;
  }

  constexpr Vector2D& operator+=(const std::array<int, 2>& Other) noexcept
  {
    x += Other[0];
    y += Other[1];
    return *this;
  }

  constexpr Vector2D& operator-=(const std::array<int, 2>& Other) noexcept
  {
    x -= Other[0];
    y -= Other[1];
    return *this;
  }

  constexpr Vector2D& operator*=(const std::array<int, 2>& Other) noexcept
  {
    x *= Other[0];
    y *= Other[1];
    return *this;
  }

  constexpr Vector2D& operator/=(const std::array<int, 2>& Other) noexcept
  {
    x /= Other[0];
    y /= Other[1];
    return *this;
  }

  constexpr Vector2D& operator+=(value_type S) noexcept
  {
    x += S;
    y += S;
    return *this;
  }

  constexpr Vector2D& operator-=(value_type S) noexcept
  {
    x -= S;
    y -= S;
    return *this;
  }

  constexpr Vector2D& operator*=(value_type S) noexcept
  {
    x *= S;
    y *= S;
    return *this;
  }

  constexpr Vector2D& operator/=(value_type S) noexcept
  {
    return *this *= (static_cast<value_type>(1.0) / S);
  }

  constexpr bool operator==(const Vector2D& V) const noexcept
  {
    return V.x == x && V.y == y;
  }

  constexpr bool operator!=(const Vector2D& V) const noexcept
  {
    return V.x != x || V.y != y;
  }

  constexpr Vector2D Clamped(const Vector2D& Min, const Vector2D& Max) const noexcept
  {
    return {x < Min.x ? Min.x : x > Max.x ? Max.x : x,
            y < Min.y ? Min.y : y > Max.y ? Max.y : y};
  }

  value_type Length() const noexcept { return std::sqrt(LengthSq()); }

  value_type LengthInv() const noexcept
  {
    return static_cast<value_type>(1.0) / Length();
  }

  constexpr value_type LengthSq() const noexcept { return Dot(*this); }

  constexpr value_type Dot(const Vector2D& Other) const noexcept
  {
    return x * Other.x + y * Other.y;
  }

  constexpr value_type Cross(const Vector2D& Other) const noexcept
  {
    return x * Other.y - y * Other.x;
  }

  constexpr Vector2D Projection(const Vector2D& Onto) const noexcept
  {
    return Onto.LengthSq() ? Onto * Dot(Onto) / Onto.LengthSq() : Zero();
  }

  Vector2D Normalized() const noexcept { return *this / Length(); }

  Vector2D& Normalize() noexcept { return *this /= Length(); }

  Vector2D RotatedAng(const value_type Angle) const noexcept
  {
    constexpr value_type r = 3.14159 / 180.0;
    const value_type radian = Angle * r;
    return {RotatedRad(radian)};
  }

  Vector2D RotatedRad(const value_type Radian) const noexcept
  {
    const value_type s = std::sin(Radian);
    const value_type c = std::cos(Radian);
    return {x * c - y * s, x * s + y * c};
  }

  Vector2D& RotateAng(const value_type Angle) noexcept
  {
    return *this = RotatedAng(Angle);
  }

  Vector2D& RotateRad(const value_type Radian) noexcept
  {
    return *this = RotatedRad(Radian);
  }

  constexpr bool IsZero() const noexcept
  {
    return x == static_cast<value_type>(0.0) && y == static_cast<value_type>(0.0);
  }

  constexpr std::array<int, 2> Int() const noexcept
  {
    return {static_cast<int>(x), static_cast<int>(y)};
  }

  static constexpr Vector2D Zero() noexcept { return {0, 0}; }

  static constexpr Vector2D One() noexcept { return {1, 1}; }

  static constexpr Vector2D UnitX() noexcept { return {1, 0}; }

  static constexpr Vector2D UnitY() noexcept { return {0, 1}; }
};

template <class Type, class U>
inline constexpr Vector2D<Type> operator+(U S, const Vector2D<Type>& V) noexcept
{
  return {static_cast<Type>(S) + V.x, static_cast<Type>(S) + V.y};
}

template <class Type, class U>
inline constexpr Vector2D<Type> operator-(U S, const Vector2D<Type>& V) noexcept
{
  return {static_cast<Type>(S) - V.x, static_cast<Type>(S) - V.y};
}

template <class Type, class U>
inline constexpr Vector2D<Type> operator*(U S, const Vector2D<Type>& V) noexcept
{
  return {static_cast<Type>(S) * V.x, static_cast<Type>(S) * V.y};
}

template <class Type, class U>
inline constexpr Vector2D<Type> operator/(U S, const Vector2D<Type>& V) noexcept
{
  return {static_cast<Type>(S) / V.x, static_cast<Type>(S) / V.y};
}

template <class CharType, class Type>
inline std::basic_ostream<CharType>& operator<<(std::basic_ostream<CharType>& OUT,
                                                const Vector2D<Type>& V)
{
  return OUT << V.x << CharType(',') << V.y;
}

using Float2 = Vector2D<float>;
using Vec2 = Vector2D<double>;

template <class Type>
struct Vector3D
{
  template <class U>
  using vector_type = Vector3D<U>;
  using value_type = Type;
  value_type x, y, z;

  Vector3D() noexcept = default;
  constexpr Vector3D(const Vector3D&) noexcept = default;

  template <class X, class Y, class Z>
  constexpr Vector3D(X _X, Y _Y, Z _Z) noexcept
    : x(static_cast<value_type>(_X)),
      y(static_cast<value_type>(_Y)),
      z(static_cast<value_type>(_Z))
  {
  }

  constexpr Vector3D(value_type _X, value_type _Y, value_type _Z) noexcept
    : x(_X), y(_Y), z(_Z)
  {
  }

  template <class U>
  constexpr Vector3D(const Vector3D<U>& V) noexcept
    : x(static_cast<value_type>(V.x)),
      y(static_cast<value_type>(V.y)),
      z(static_cast<value_type>(V.z))
  {
  }

  constexpr Vector3D& Set(value_type _X, value_type _Y, value_type _Z) noexcept
  {
    x = _X;
    y = _Y;
    z = _Z;
    return *this;
  }

  constexpr Vector3D& Set(const Vector3D& V) noexcept { return *this = V; }

  constexpr Vector3D operator+() const noexcept { return *this; }

  constexpr Vector3D operator-() const noexcept { return {-x, -y, -z}; }

  constexpr Vector3D operator+(const Vector3D& Other) const noexcept
  {
    return {x + Other.x, y + Other.y, z + Other.z};
  }

  constexpr Vector3D operator-(const Vector3D& Other) const noexcept
  {
    return {x - Other.x, y - Other.y, z - Other.z};
  }

  constexpr Vector3D operator*(const Vector3D& Other) const noexcept
  {
    return {x * Other.x, y * Other.y, z * Other.z};
  }

  constexpr Vector3D operator/(const Vector3D& Other) const noexcept
  {
    return {x / Other.x, y / Other.y, z / Other.z};
  }

  constexpr Vector3D operator+(const std::array<int, 3>& Other) const noexcept
  {
    return {x + Other[0], y + Other[1], z + Other[2]};
  }

  constexpr Vector3D operator-(const std::array<int, 3>& Other) const noexcept
  {
    return {x - Other[0], y - Other[1], z - Other[2]};
  }

  constexpr Vector3D operator*(const std::array<int, 3>& Other) const noexcept
  {
    return {x * Other[0], y * Other[1], z * Other[2]};
  }

  constexpr Vector3D operator/(const std::array<int, 3>& Other) const noexcept
  {
    return {x / Other[0], y / Other[1], z / Other[2]};
  }

  constexpr Vector3D operator+(value_type S) const noexcept
  {
    return {x + S, y + S, z + S};
  }

  constexpr Vector3D operator-(value_type S) const noexcept
  {
    return {x - S, y - S, z - S};
  }

  constexpr Vector3D operator*(value_type S) const noexcept
  {
    return {x * S, y * S, z * S};
  }

  constexpr Vector3D operator/(value_type S) const noexcept
  {
    return *this * (static_cast<value_type>(1.0) / S);
  }

  constexpr Vector3D& operator+=(const Vector3D& Other) noexcept
  {
    x += Other.x;
    y += Other.y;
    z += Other.z;
    return *this;
  }

  constexpr Vector3D& operator-=(const Vector3D& Other) noexcept
  {
    x -= Other.x;
    y -= Other.y;
    z -= Other.z;
    return *this;
  }

  constexpr Vector3D& operator*=(const Vector3D& Other) noexcept
  {
    x *= Other.x;
    y *= Other.y;
    z *= Other.z;
    return *this;
  }

  constexpr Vector3D& operator/=(const Vector3D& Other) noexcept
  {
    x /= Other.x;
    y /= Other.y;
    z /= Other.z;
    return *this;
  }

  constexpr Vector3D& operator+=(const std::array<int, 3>& Other) noexcept
  {
    x += Other[0];
    y += Other[1];
    z += Other[2];
    return *this;
  }

  constexpr Vector3D& operator-=(const std::array<int, 3>& Other) noexcept
  {
    x -= Other[0];
    y -= Other[1];
    z -= Other[2];
    return *this;
  }

  constexpr Vector3D& operator*=(const std::array<int, 3>& Other) noexcept
  {
    x *= Other[0];
    y *= Other[1];
    z *= Other[2];
    return *this;
  }

  constexpr Vector3D& operator/=(const std::array<int, 3>& Other) noexcept
  {
    x /= Other[0];
    y /= Other[1];
    z /= Other[2];
    return *this;
  }

  constexpr Vector3D& operator+=(value_type S) noexcept
  {
    x += S;
    y += S;
    z += S;
    return *this;
  }

  constexpr Vector3D& operator-=(value_type S) noexcept
  {
    x -= S;
    y -= S;
    z -= S;
    return *this;
  }

  constexpr Vector3D& operator*=(value_type S) noexcept
  {
    x *= S;
    y *= S;
    z *= S;
    return *this;
  }

  constexpr Vector3D& operator/=(value_type S) noexcept
  {
    return *this *= (static_cast<value_type>(1.0) / S);
  }

  constexpr bool operator==(const Vector3D& V) const noexcept
  {
    return V.x == x && V.y == y && V.z == z;
  }

  constexpr bool operator!=(const Vector3D& V) const noexcept
  {
    return V.x != x || V.y != y || V.z != z;
  }

  constexpr Vector3D Clamped(const Vector3D& Min, const Vector3D& Max) const noexcept
  {
    return {x < Min.x ? Min.x : x > Max.x ? Max.x : x,
            y < Min.y ? Min.y : y > Max.y ? Max.y : y,
            z < Min.z ? Min.z : z > Max.z ? Max.z : z};
  }

  Vector3D RotatedAng(const Vector3D& Unit_Vec3, const value_type Ang) const noexcept
  {
    constexpr value_type r = 3.14159 / 180.0;
    const value_type rad = Ang * r;
    return {RotatedRad(Unit_Vec3, rad)};
  }

  Vector3D RotatedRad(const Vector3D& Unit_Vec3, const value_type Rad) const noexcept
  {
    const value_type s = std::sin(Rad);
    const value_type c = std::cos(Rad);
    return {((*this) * c) + (Unit_Vec3.Cross(*this) * s) +
            ((Unit_Vec3 * (Unit_Vec3.Dot(*this))) * (1 - c))};
  }

  Vector3D& RotateAng(const Vector3D& Unit_Vec3, const value_type Ang) noexcept
  {
    return *this = RotatedAng(Unit_Vec3, Ang);
  }

  Vector3D& RotateRad(const Vector3D& Unit_Vec3, const value_type Rad) noexcept
  {
    return *this = RotatedRad(Unit_Vec3, Rad);
  }

  constexpr bool IsZero() const noexcept
  {
    return x == static_cast<value_type>(0.0) && y == static_cast<value_type>(0.0) &&
           z == static_cast<value_type>(0.0);
  }

  constexpr value_type Dot(const Vector3D& Other) const noexcept
  {
    return x * Other.x + y * Other.y + z * Other.z;
  }

  constexpr Vector3D Cross(const Vector3D& Other) const noexcept
  {
    return {y * Other.z - z * Other.y, z * Other.x - x * Other.z,
            x * Other.y - y * Other.x};
  }

  value_type Length() const noexcept { return std::sqrt(LengthSq()); }

  value_type LengthInv() const noexcept
  {
    return static_cast<value_type>(1.0) / Length();
  }

  constexpr value_type LengthSq() const noexcept { return Dot(*this); }

  Vector3D Normalized() const { return *this / Length(); }

  Vector3D& Normalize() { return *this /= Length(); }

  constexpr Vector2D<value_type> XY() const noexcept { return {x, y}; }

  constexpr Vector2D<value_type> XZ() const noexcept { return {x, z}; }

  constexpr Vector2D<value_type> YZ() const noexcept { return {y, z}; }

  constexpr std::array<int, 3> Int() const noexcept
  {
    return {static_cast<int>(x), static_cast<int>(y), static_cast<int>(z)};
  }

  constexpr std::array<int, 2> IntXY() const noexcept
  {
    return {static_cast<int>(x), static_cast<int>(y)};
  }

  constexpr std::array<int, 2> IntXZ() const noexcept
  {
    return {static_cast<int>(x), static_cast<int>(z)};
  }

  constexpr std::array<int, 2> IntYZ() const noexcept
  {
    return {static_cast<int>(y), static_cast<int>(z)};
  }

  static constexpr Vector3D Zero() { return {0, 0, 0}; }

  static constexpr Vector3D One() { return {1, 1, 1}; }

  static constexpr Vector3D UnitX() { return {1, 0, 0}; }

  static constexpr Vector3D UnitY() { return {0, 1, 0}; }

  static constexpr Vector3D UnitZ() { return {0, 0, 1}; }
};

template <class Type, class U>
inline constexpr Vector3D<Type> operator+(U S, const Vector3D<Type>& V) noexcept
{
  return {static_cast<Type>(S) + V.x, static_cast<Type>(S) + V.y,
          static_cast<Type>(S) + V.z};
}

template <class Type, class U>
inline constexpr Vector3D<Type> operator-(U S, const Vector3D<Type>& V) noexcept
{
  return {static_cast<Type>(S) - V.x, static_cast<Type>(S) - V.y,
          static_cast<Type>(S) - V.z};
}

template <class Type, class U>
inline constexpr Vector3D<Type> operator*(U S, const Vector3D<Type>& V) noexcept
{
  return {static_cast<Type>(S) * V.x, static_cast<Type>(S) * V.y,
          static_cast<Type>(S) * V.z};
}

template <class Type, class U>
inline constexpr Vector3D<Type> operator/(U S, const Vector3D<Type>& V) noexcept
{
  return {static_cast<Type>(S) / V.x, static_cast<Type>(S) / V.y,
          static_cast<Type>(S) / V.z};
}

template <class CharType, class Type>
inline std::basic_ostream<CharType>& operator<<(std::basic_ostream<CharType>& OUT,
                                                const Vector3D<Type>& V)
{
  return OUT << V.x << CharType(',') << V.y << CharType(',') << V.z;
}

using Float3 = Vector3D<float>;
using Vec3 = Vector3D<double>;
}  // namespace dem
