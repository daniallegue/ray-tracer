#include "draw.h"
#include "common.h"
#include <framework/opengl_includes.h>
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#ifdef __APPLE__
#include <OpenGL/GLU.h>
#else
#ifdef WIN32
// Windows.h includes a ton of stuff we don't need, this macro tells it to include less junk.
#define WIN32_LEAN_AND_MEAN
// Disable legacy macro of min/max which breaks completely valid C++ code (std::min/std::max won't work).
#define NOMINMAX
// GLU requires Windows.h on Windows :-(.
#include <Windows.h>
#endif
#include <GL/glu.h>
#endif
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/mat4x4.hpp>
DISABLE_WARNINGS_POP()
#include <algorithm>

bool enableDebugDraw = false;
bool enableBlurDebug; 

static void setMaterial(const Material& material)
{
    // Set the material color of the shape.
    const glm::vec4 kd4 { material.kd, 1.0f };
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, glm::value_ptr(kd4));

    const glm::vec4 zero { 0.0f };
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, glm::value_ptr(zero));
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, glm::value_ptr(zero));
}

void drawLine(const glm::vec3 origin, const glm::vec3 line, const glm::vec3& color)
{
    glBegin(GL_LINES);
    glColor3fv(glm::value_ptr(color));
    glVertex3fv(glm::value_ptr(origin));
    glm::vec3 end = { origin.x + line.x, origin.y + line.y, origin.z + origin.z };
    glVertex3fv(glm::value_ptr(end));
    glEnd();
}

void drawSegment(const glm::vec3 origin, const glm::vec3 end, const glm::vec3& color)
{
    glBegin(GL_LINES);
    glColor3fv(glm::value_ptr(color));
    glVertex3fv(glm::value_ptr(origin));
    glVertex3fv(glm::value_ptr(end));
    glLineWidth(1000.0f);
    glEnd();
}

void drawExampleOfCustomVisualDebug()
{
    glBegin(GL_TRIANGLES);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    glEnd();
}


void drawTriangle (const Vertex& v0, const Vertex& v1, const Vertex& v2 ) {
    glBegin(GL_TRIANGLES);
        glNormal3fv(glm::value_ptr(v0.normal));
        glVertex3fv(glm::value_ptr(v0.position));
        glNormal3fv(glm::value_ptr(v1.normal));
        glVertex3fv(glm::value_ptr(v1.position));
        glNormal3fv(glm::value_ptr(v2.normal));
        glVertex3fv(glm::value_ptr(v2.position));
    glEnd();
}

void drawMesh(const Mesh& mesh)
{
    setMaterial(mesh.material);

    glBegin(GL_TRIANGLES);
    for (const auto& triangleIndex : mesh.triangles) {
        for (int i = 0; i < 3; i++) {
            const auto& vertex = mesh.vertices[triangleIndex[i]];
            glNormal3fv(glm::value_ptr(vertex.normal));
            glVertex3fv(glm::value_ptr(vertex.position));
        }
    }
    glEnd();
}

static void drawSphereInternal(const glm::vec3& center, float radius)
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    const glm::mat4 transform = glm::translate(glm::identity<glm::mat4>(), center);
    glMultMatrixf(glm::value_ptr(transform));
    auto quadric = gluNewQuadric();
    gluSphere(quadric, radius, 40, 20);
    gluDeleteQuadric(quadric);
    glPopMatrix();
}

void drawBezierCurve(glm::vec3 pos, float radius, glm::vec3 color)
{
    //Use same points as other function 
    //Draw first vertex
    //glVertex3f(pos.x, pos.y, pos.z);
    int counter = 0;
    for (float i = 0; i < 1.3f; i += 0.01f) {
        // Generate deterministic samples
        float u = 1.0f - i;
        float u2 = glm::pow(u, 2);
        float u3 = glm::pow(u, 3);
        float t2 = glm::pow(i, 2);
        float t3 = glm::pow(i, 3);
        glm::vec3 p0 = (glm::vec3(0, 0, 0) * 1.0f) + pos; 
        glm::vec3 p1 = (glm::vec3(1, 2, 2) * 1.0f) + pos;
        glm::vec3 p2 = (glm::vec3(1, 2, 2) * 1.0f) + pos;
        glm::vec3 p3 = (glm::vec3(3, 1, 0) * 1.0f) + pos;
        glm::vec3 newPos = (u3 * p0) + (3.0f * u2 * i * p1) + (3.0f * u * t2 * p2) + (t3 * p3);
        //glVertex3f(newPos.x, newPos.y, newPos.z);
        //Generate a sphere every 10 samples
        if (counter % 10 == 0) {
            drawSphere(newPos, radius, color, 0.3f);
        }
        counter++;
    }

}

void drawSphere(const Sphere& sphere)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    setMaterial(sphere.material);
    drawSphereInternal(sphere.center, sphere.radius);
    glPopAttrib();
}

void drawSphere(const glm::vec3& center, float radius, const glm::vec3& color, float opacity)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glColor4f(color.r, color.g, color.b, opacity);
    glPolygonMode(GL_FRONT, GL_FILL);
    glPolygonMode(GL_BACK, GL_FILL);
    drawSphereInternal(center, radius);
    glPopAttrib();
}

void drawSphere(const glm::vec3& center, float radius, const glm::vec3& color /*= glm::vec3(1.0f)*/)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glColor4f(color.r, color.g, color.b, 1.0f);
    glPolygonMode(GL_FRONT, GL_FILL);
    glPolygonMode(GL_BACK, GL_FILL);
    drawSphereInternal(center, radius);
    glPopAttrib();
}

static void drawAABBInternal(const AxisAlignedBox& box)
{
    glPushMatrix();

    // front      back
    // 3 ----- 2  7 ----- 6
    // |       |  |       |
    // |       |  |       |
    // 0 ------1  4 ------5

    glBegin(GL_QUADS);
    glNormal3f(0, 0, -1);
    glVertex3f(box.lower.x, box.upper.y, box.lower.z); //3
    glVertex3f(box.upper.x, box.upper.y, box.lower.z); //2
    glVertex3f(box.upper.x, box.lower.y, box.lower.z); //1
    glVertex3f(box.lower.x, box.lower.y, box.lower.z); //0

    glNormal3f(0, 0, 1);
    glVertex3f(box.upper.x, box.lower.y, box.upper.z); //5
    glVertex3f(box.upper.x, box.upper.y, box.upper.z); //6
    glVertex3f(box.lower.x, box.upper.y, box.upper.z); //7
    glVertex3f(box.lower.x, box.lower.y, box.upper.z); //4

    glNormal3f(1, 0, 0);
    glVertex3f(box.upper.x, box.upper.y, box.lower.z); //2
    glVertex3f(box.upper.x, box.upper.y, box.upper.z); //6
    glVertex3f(box.upper.x, box.lower.y, box.upper.z); //5
    glVertex3f(box.upper.x, box.lower.y, box.lower.z); //1

    glNormal3f(-1, 0, 0);
    glVertex3f(box.lower.x, box.lower.y, box.upper.z); //4
    glVertex3f(box.lower.x, box.upper.y, box.upper.z); //7
    glVertex3f(box.lower.x, box.upper.y, box.lower.z); //3
    glVertex3f(box.lower.x, box.lower.y, box.lower.z); //0

    glNormal3f(0, 1, 0);
    glVertex3f(box.lower.x, box.upper.y, box.upper.z); //7
    glVertex3f(box.upper.x, box.upper.y, box.upper.z); //6
    glVertex3f(box.upper.x, box.upper.y, box.lower.z); //2
    glVertex3f(box.lower.x, box.upper.y, box.lower.z); //3

    glNormal3f(0, -1, 0);
    glVertex3f(box.upper.x, box.lower.y, box.lower.z); //1
    glVertex3f(box.upper.x, box.lower.y, box.upper.z); //5
    glVertex3f(box.lower.x, box.lower.y, box.upper.z); //4
    glVertex3f(box.lower.x, box.lower.y, box.lower.z); //0
    glEnd();

    glPopMatrix();
}

void drawAABB(const AxisAlignedBox& box, DrawMode drawMode, const glm::vec3& color, float transparency)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glColor4f(color.r, color.g, color.b, transparency);
    if (drawMode == DrawMode::Filled) {
        glPolygonMode(GL_FRONT, GL_FILL);
        glPolygonMode(GL_BACK, GL_FILL);
    } else {
        glPolygonMode(GL_FRONT, GL_LINE);
        glPolygonMode(GL_BACK, GL_LINE);
    }
    drawAABBInternal(box);
    glPopAttrib();
}

void drawScene(const Scene& scene)
{
    for (const auto& mesh : scene.meshes)
        drawMesh(mesh);
    for (const auto& sphere : scene.spheres)
        drawSphere(sphere);
}

void drawRay(const Ray& ray, const glm::vec3& color)
{
    if (!enableDebugDraw)
        return;

    const glm::vec3 hitPoint = ray.origin + std::clamp(ray.t, 0.0f, 100.0f) * ray.direction;
    const bool hit = (ray.t < std::numeric_limits<float>::max());

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);

    glColor3fv(glm::value_ptr(color));
    glVertex3fv(glm::value_ptr(ray.origin));
    glColor3fv(glm::value_ptr(color));
    glVertex3fv(glm::value_ptr(hitPoint));
    glEnd();

    if (hit)
        drawSphere(hitPoint, 0.005f, color);

    glPopAttrib();
}
