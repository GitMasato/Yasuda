# pragma once
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

namespace dem
{

	struct DataIO
	{

		template<typename T>
		void TransformString(const std::string& S, T& V) const noexcept
		{
			T t;
			std::stringstream ss;
			ss << S;
			ss >> t;
			V = t;
		}


		int StringToInt(const std::string& Str) const noexcept
		{
			int t;
			std::stringstream ss;
			ss << Str;
			ss >> t;
			return t;
		}


		double StringToDouble(const std::string& Str) const noexcept
		{
			double t;
			std::stringstream ss;
			ss << Str;
			ss >> t;
			return t;
		}


		bool StringToBool(const std::string& Str) const noexcept
		{
			bool t;
			std::stringstream ss;
			ss << Str;
			ss >> t;
			return t;
		}


		void ErrorInput(const std::string& Error_file) const noexcept
		{
			std::cerr << "Error to read [" << Error_file << "]" << std::endl << std::endl;
			exit(1);
		}


		void ErrorOutput(const std::string& Error_file) const noexcept
		{
			std::cerr << "Error to output [" << Error_file << "]" << std::endl << std::endl;
			exit(1);
		}


		void SplitLine(const std::string& Line, std::vector<std::string>& Split_Line) const noexcept
		{
			std::string v;
			std::stringstream ss(Line);

			while (ss >> v)
			{
				const std::size_t num = v.find('#');
				if (num != std::string::npos) v.erase(num);
				Split_Line.push_back(v);
			}
		}


		std::list<std::string> Split(const std::string& Line, const std::string& Delim) const noexcept
		{
			std::list<std::string> result;
			std::string Str = Line;
			std::size_t cutAt;

			while ((cutAt = static_cast<std::size_t>(Str.find_first_of(Delim))) != Str.npos)
			{
				if (cutAt > 0) { result.push_back(Str.substr(0, cutAt)); }
				const std::size_t substr_temp = cutAt + 1;
				Str = Str.substr(substr_temp);
			}

			if (Str.length() > 0) { result.push_back(Str); }
			return result;
		}


	};

}