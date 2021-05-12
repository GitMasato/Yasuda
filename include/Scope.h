#pragma once
#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "Vector.h"

namespace dem
{
struct Scope
{
  template <typename T>
  void ShowVariable(const T V) const noexcept
  {
    std::cout << "(" << V << ")" << std::endl;
  }

  template <typename T>
  void ShowVariable(const T V, const std::string& Message) const noexcept
  {
    std::cout << Message << " (" << V << ")" << std::endl;
  }

  template <typename T>
  void ShowArray(const T& A) const noexcept
  {
    for (const auto& v : A)
    {
      std::size_t i = &v - &A[0];
      std::cout << "(" << i << "," << v << ")" << std::endl;
    }
  }

  template <typename T>
  void ShowArray(const T& A, const std::string& Message) const noexcept
  {
    for (const auto& v : A)
    {
      std::size_t i = &v - &A[0];
      std::cout << Message << " (" << i << "," << v << ")" << std::endl;
    }
  }

  template <typename T>
  void ShowMultipleArray(const T& M) const noexcept
  {
    for (const auto& w : M)
    {
      std::size_t k = &w - &M[0];

      for (const auto& v : w)
      {
        std::size_t i = &v - &w[0];
        std::cout << "(" << k << "," << i << "," << v << ")" << std::endl;
      }
    }
  }

  template <typename T>
  void ShowMultipleArray(const T& M, const std::string& Message) const noexcept
  {
    for (const auto& w : M)
    {
      std::size_t k = &w - &M[0];

      for (const auto& v : w)
      {
        std::size_t i = &v - &w[0];
        std::cout << Message << " (" << k << "," << i << "," << v << ")"
                  << std::endl;
      }
    }
  }
};
}  // namespace dem
