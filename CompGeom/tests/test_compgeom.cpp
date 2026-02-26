#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include "Geometry.hpp"

void test_signed_area() {
    std::cout << "Testing JGeometry::getSignedArea..." << std::endl;
    
    // A unit square: (0,0), (1,0), (1,1), (0,1)
    double x[] = {0.0, 1.0, 1.0, 0.0};
    double y[] = {0.0, 0.0, 1.0, 1.0};
    
    double area = JGeometry::getSignedArea(x, y, 4);
    
    // Area should be 1.0 (CCW) or -1.0 (CW)
    assert(std::abs(std::abs(area) - 1.0) < 1e-9);
    
    std::cout << "JGeometry::getSignedArea passed!" << std::endl;
}

int main() {
    test_signed_area();
    std::cout << "All CompGeom tests passed!" << std::endl;
    return 0;
}
