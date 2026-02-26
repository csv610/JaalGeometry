#include <iostream>
#include <vector>
#include "Geometry.hpp"

/**
 * Example: Point-in-Polygon Test using JGeometry
 */
int main() {
    // 1. Define a simple square polygon (0,0) to (1,1)
    std::vector<Point2D> polygon;
    polygon.push_back({0.0, 0.0});
    polygon.push_back({1.0, 0.0});
    polygon.push_back({1.0, 1.0});
    polygon.push_back({0.0, 1.0});

    // 2. Define query points
    Point2D insidePoint = {0.5, 0.5};
    Point2D outsidePoint = {1.5, 1.5};

    // 3. Perform the tests
    int result1 = JGeometry::isInside(polygon, insidePoint);
    int result2 = JGeometry::isInside(polygon, outsidePoint);

    // 4. Print results
    std::cout << "Point (0.5, 0.5) is " << (result1 ? "INSIDE" : "OUTSIDE") << " the polygon." << std::endl;
    std::cout << "Point (1.5, 1.5) is " << (result2 ? "INSIDE" : "OUTSIDE") << " the polygon." << std::endl;

    // 5. Calculate area
    double area = JGeometry::getSignedArea(polygon);
    std::cout << "Polygon Area: " << area << std::endl;

    return 0;
}
