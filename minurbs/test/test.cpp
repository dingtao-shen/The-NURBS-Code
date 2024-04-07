#include<iostream>
#include <vector>
#include<algorithm>
#include <cmath>
#include <limits>
#include <cstddef>
#include <iomanip>
#include <type_traits>
using namespace std;

int main(){
    cout.precision(16);
    
    double eps = std::numeric_limits<double>::epsilon();
    cout << eps << endl;

	return 0;
}