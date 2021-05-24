# pragma once
#include <iostream>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <list>
#include <algorithm>

using namespace std;

class DataIO
{
private:

public:

	int int_string(const string& str){				//string型をint型へ変換する関数			

		int t;
		stringstream ss;
		ss << str;
		ss >> t;
		return t;
	}

	double double_string(const string& str){		//string型をdouble型へ変換する関数

		double t;
		stringstream ss;
		ss << str;
		ss >> t;
		return t;
	}

	bool bool_string(const string& str){			//string型をbool型へ変換する関数

		bool t;
		stringstream ss;
		ss << str;
		ss >> t;
		return t;
	}

	list<string> split(string str, string delim){	//split関数

		list<string> result;
		int cutAt;

		while ((cutAt = (int)(str.find_first_of(delim))) != str.npos){
			if (cutAt > 0){
				result.push_back(str.substr(0, cutAt));
			}
			str = str.substr(cutAt + 1);
		}
		if (str.length() > 0){

			result.push_back(str);
		}
		return result;
	}

};
