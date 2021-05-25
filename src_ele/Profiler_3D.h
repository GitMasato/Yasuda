# pragma once
#include <iostream>
#include <chrono>

class Profiler
{
private:

	decltype(std::chrono::high_resolution_clock::now()) begin
		= std::chrono::high_resolution_clock::now();

public:

	long long elapsed() const
	{
		const auto end = std::chrono::high_resolution_clock::now();

		return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	}

	int print() const
	{
		return (int)(0.001 * elapsed());
	}
};
