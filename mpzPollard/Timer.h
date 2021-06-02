#pragma once
#include <chrono>

class Timer {
	std::chrono::steady_clock::time_point begintime;
	std::chrono::steady_clock::time_point endtime;
public:
	void start() {
		begintime = std::chrono::steady_clock::now();
	};
	void end() {
		endtime = std::chrono::steady_clock::now();
	}
	int64_t ms() {
		return std::chrono::duration_cast<std::chrono::milliseconds>(endtime - begintime).count();
	}
	int64_t microsec() {
		return std::chrono::duration_cast<std::chrono::microseconds>(endtime - begintime).count();
	}
};