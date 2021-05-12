#pragma once
#include <chrono>
#include <ctime>
#include <sstream>
#include <string>

namespace dem
{
class Clock
{
private:
  decltype(std::chrono::high_resolution_clock::now()) begin =
    std::chrono::high_resolution_clock::now();

public:
  long long Elapsed() const noexcept
  {
    const auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
      .count();
  }

  int GetElapsedTime() const noexcept
  {
    return static_cast<int>(0.001 * Elapsed());
  }

  std::string GetLocalDate() const noexcept
  {
    time_t now = time(nullptr);
    std::string date = ctime(&now);
    return DeleteNewLine(date);
  }

  std::string DeleteNewLine(const std::string& Target) const noexcept
  {
    constexpr char CR = '\r';
    constexpr char LF = '\n';
    std::string str;

    for (const auto& s : Target)
    {
      {
        if ((s != CR) && (s != LF))
        {
          str += s;
        }
      }
    }
    return str;
  }
};

}  // namespace dem
