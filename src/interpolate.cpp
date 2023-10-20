#include "interpolate.h"
#include <glm/geometric.hpp>

// TODO Standard feature
// Given three triangle vertices and a point on the triangle, compute the corresponding barycentric coordinates of the point.
// and return a vec3 with the barycentric coordinates (alpha, beta, gamma).
// - v0;     Triangle vertex 0
// - v1;     Triangle vertex 1
// - v2;     Triangle vertex 2
// - p;      Point on triangle
// - return; Corresponding barycentric coordinates for point p.
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    glm::vec3 a = v1 - v0;
    glm::vec3 b = v2 - v0;
    glm::vec3 c = p - v0;

    float beta = (glm::dot(a, a) * glm::dot(a, c) - glm::dot(a, b) * glm::dot(b, c)) / (glm::dot(a, a) * glm::dot(b, b) - glm::dot(a, b) * glm::dot(a, b));
    float gamma = (glm::dot(a, a) * glm::dot(b, c) - glm::dot(a, b) * glm::dot(a, c)) / (glm::dot(a, a) * glm::dot(b, b) - glm::dot(a, b) * glm::dot(a, b));
    float alpha = 1.0f - beta - gamma;

    return glm::vec3 { alpha, beta, gamma };
}

// TODO Standard feature
// Linearly interpolate three normals using barycentric coordinates.
// - n0;     Triangle normal 0
// - n1;     Triangle normal 1
// - n2;     Triangle normal 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated normal.
// This method is unit-tested, so do not change the function signature.
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 bc)
{
    glm::vec3 interpolated = bc.x * n0 + bc.y * n1 + bc.z * n2;
    return interpolated;

}

// TODO Standard feature
// Linearly interpolate three texture coordinates using barycentric coordinates.
// - n0;     Triangle texture coordinate 0
// - n1;     Triangle texture coordinate 1
// - n2;     Triangle texture coordinate 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated texturre coordinate.
// This method is unit-tested, so do not change the function signature.
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 bc)
{
    glm::vec2 interpolated = bc.x * t0 + bc.y * t1 + bc.z * t2;
    return interpolated;
}
