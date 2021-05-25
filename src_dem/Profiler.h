# pragma once
#include <chrono>

namespace dem
{
	class Profiler
	{
	private:

		decltype(std::chrono::high_resolution_clock::now()) begin
			= std::chrono::high_resolution_clock::now();

	public:

		long long Elapsed() const noexcept
		{
			const auto end = std::chrono::high_resolution_clock::now();

			return std::chrono::
				duration_cast<std::chrono::milliseconds>(end - begin).count();
		}

		int Print() const noexcept
		{
			return static_cast<int>(0.001 * Elapsed());
		}

	};

}
