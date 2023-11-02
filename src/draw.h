#pragma once
#include "scene.h"
#include <framework/mesh.h>
#include <framework/ray.h>
#include <utility> // std::forward

// Flag to enable/disable the debug drawing.
extern bool enableDebugDraw;

// Add your own custom visual debug draw functions here then implement it in draw.cpp.
// You are free to modify the example one however you like.
//
// Please add following lines on top inside your custom draw function:
//
// if (!enableDebugDraw)
//     return;
//
void drawExampleOfCustomVisualDebug();

void drawRay(const Ray& ray, const glm::vec3& color = glm::vec3(1.0f));

void drawAABB(const AxisAlignedBox& box, DrawMode drawMode = DrawMode::Filled, const glm::vec3& color = glm::vec3(1.0f), float transparency = 1.0f);
void drawLine(const glm::vec3 origin, const glm::vec3 line, const glm::vec3& color);
void drawBezierCurveSpheres(glm::vec3 pos, float radius, glm::vec3 color);
void drawBezierCurveMeshes(Mesh mesh);
void drawSphere(const glm::vec3& center, float radius, const glm::vec3& color, float opacity);
void drawSegment(const glm::vec3 origin, const glm::vec3 end, const glm::vec3& color);
void drawTriangle (const Vertex& v0, const Vertex& v1, const Vertex& v2 );
void drawMesh(const Mesh& mesh);
void drawSphere(const Sphere& sphere);
void drawSphere(const glm::vec3& center, float radius, const glm::vec3& color = glm::vec3(1.0f));
void drawScene(const Scene& scene);


/// <summary>
/// Draws a rectangle inside an <see cref="AxisAlignedBox"/>
/// </summary>
/// <param name="aabb"></param>
/// <param name="p"></param>
/// <param name="axis"></param>
/// <param name="color"></param>
void drawSplitPlane(const AxisAlignedBox& aabb, const glm::vec3& p, const glm::length_t axis, const glm::vec4& color);

