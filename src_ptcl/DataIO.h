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

	int int_string(const string& str){				//string�^��int�^�֕ϊ�����֐�			

		int t;
		stringstream ss;
		ss << str;
		ss >> t;
		return t;
	}

	double double_string(const string& str){		//string�^��double�^�֕ϊ�����֐�

		double t;
		stringstream ss;
		ss << str;
		ss >> t;
		return t;
	}

	bool bool_string(const string& str){			//string�^��bool�^�֕ϊ�����֐�

		bool t;
		stringstream ss;
		ss << str;
		ss >> t;
		return t;
	}

	list<string> split(string str, string delim){	//split�֐�

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
