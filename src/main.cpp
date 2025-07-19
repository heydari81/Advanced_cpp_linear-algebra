#include "algebra.h"

#include <gtest/gtest.h>
#include <iostream>
using namespace std;
template<typename T>
using MATRIX = std::vector<std::vector<T>>;
int main(int argc, char **argv) {
	if (false) // make false to run unit-tests
	{
        MATRIX<int> mat = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    	return 0;
		// debug section
	} else {
		::testing::InitGoogleTest(&argc, argv);
		std::cout << "RUNNING TESTS ..." << std::endl;
		int ret{RUN_ALL_TESTS()};
		if (!ret)
			std::cout << "<<<SUCCESS>>>" << std::endl;
		else
			std::cerr << "FAILED" << std::endl;
	}
	return 0;
}
