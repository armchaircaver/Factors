#include <vector>
#include <stdio.h>
#include <stdint.h>
#include <ctime>	
#include <algorithm>  // sort
#include <string>
#include <iostream>
#include <corecrt_io.h>
#include <io.h>
#include <fcntl.h>

#include "../FactorsA/FactorsA.h"
#include "../mpzPollard/Timer.h"




void printvec(std::vector<uint64_t> v) {
	for (int i = 0; i < v.size() - 1; i++)
		std::wcout << v[i] << " * ";
	std::wcout << v[v.size() - 1] << "\n";
}

bool isNumber(const std::string& str){
	for (char const& c : str) 
		if (std::isdigit(c) == 0) return false;
	return true;
}

int main(int argc, char** argv) {

	int mode = _setmode(_fileno(stdout), _O_U16TEXT);

    // string representation of maximum allowed 64 bit number
	std::string maxintstr = std::to_string(UINT64_MAX);

	Timer tim;
	while (1) {

		std::string numberstr;
		std::wcout << "Enter number to factorise: ";
		std::cin >> numberstr;
		if (numberstr.size() == 0)
			return 0;

		if (! isNumber(numberstr) )		{
			wprintf(L"Not a number\n"); 
			continue;
		}
		
		if (numberstr.size() > 20) {
			wprintf(L"too large\n");
			continue;
		}
		
		if (numberstr.size() == 20 && numberstr.compare(maxintstr) > 0 ) {
			wprintf(L"too large\n");
			continue;
		}

		char* end;
		uint64_t n = strtoull(numberstr.c_str(), &end, 10);

		if (n < 2) {
			wprintf(L"too small\n");
			continue;
		}


		tim.start();
		std::vector<uint64_t> v = factorise(n);
		tim.end();

		if (v.size() == 1)
			std::wcout << "\n" << n << " is prime \n";
		else {
			std::sort(v.begin(), v.end());
			std::wcout << "\n"<< n << " = ";
			printvec(v);
		}
		std::wcout << "\n" << tim.microsec() <<  L" μs\n\n";
	}
}