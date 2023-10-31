#include "extra.h"
#include "bvh.h"
#include "intersect.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include "texture.h"
#include <framework/trackball.h>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

//Define function for spline transformation with a cubic Bézier curve
glm::mat4 bezierTransformation(const Features& features, glm::vec3 pos, float time)
{
    float u = 1.0f - time;
    float u2 = glm::pow(u, 2);
    float u3 = glm::pow(u, 3);
    float t2 = glm::pow(time, 2);
    float t3 = glm::pow(time, 3);


    float movement = features.extra.movementFactor;
    glm::vec3 p0 = (glm::vec3(0, 0, 0) * movement) + pos;
    glm::vec3 p1 = (glm::vec3(1, 2, 2) * movement) + pos;
    glm::vec3 p2 = (glm::vec3(1, 2, 2) * movement) + pos;
    glm::vec3 p3 = (glm::vec3(3, 1, 0) * movement) + pos;
    glm::vec3 newPos = (u3 * p0) + (3.0f * u2 * time * p1) + (3.0f * u * t2 * p2) + (t3 * p3);

    //Translate identity matrix by Bezier transformation
    return glm::translate(glm::mat4(1.0f), newPos);
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

    for (int y = 0; y < screen.resolution()[1]; y++) {
        for (int x = 0; x != screen.resolution()[0]; x++) {

            glm::vec3 sceneLight; 
            RenderState state = {
                .scene = scene,
                .features = features,
                .bvh = bvh,
                .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
            };

            for (int i = 0; i < features.extra.numBlurSamples; i++) {
                //Generate sample for iteration
                float time = state.sampler.next_1d();

                //Apply to meshes
                std::vector<Mesh> meshes;
                for (int i = 0; i < scene.meshes.size(); i++) {
                    Mesh mesh = scene.meshes[i];
                    std::vector<Vertex> vertices = mesh.vertices;
                    std::vector<Vertex> updatedVertices;
                    for (int u = 0; u < vertices.size(); u++) {
                        Vertex v = vertices[u];
                        glm::vec3 pos = v.position;
                        glm::mat4 transform = bezierTransformation(features, pos, time);
                        glm::vec4 newPos = transform * glm::vec4 { pos.x, pos.y, pos.z, 1 };
                        glm::vec3 posScaled = { newPos[0] / newPos[3], newPos[1] / newPos[3], newPos[2] / newPos[3] };

                        Vertex updatedVertex = {
                            .position = posScaled,
                            .normal = v.normal,
                            .texCoord = v.texCoord
                        };
                        updatedVertices.push_back(updatedVertex);
                    }

                    Mesh newMesh = {
                        .vertices = updatedVertices,
                        .triangles = mesh.triangles,
                        .material = mesh.material
                    };
                    meshes.push_back(newMesh);
                }

                //Apply to spheres
                std::vector<Sphere> spheres;
                for (int i = 0; i < scene.spheres.size(); i++) {
                    Sphere sphere = scene.spheres[i];
                    glm::vec3 pos = sphere.center;
                    glm::mat4 transform = bezierTransformation(features, pos, time);
                    glm::vec4 newPos = transform * glm::vec4 { pos.x, pos.y, pos.z, 1 };
                    glm::vec3 posScaled = { newPos[0] / newPos[3], newPos[1] / newPos[3], newPos[2] / newPos[3] };
                    Sphere updatedSphere = {
                        .center = posScaled,
                        .radius = sphere.radius,
                        .material = sphere.material
                    };

                    spheres.push_back(updatedSphere);
                }

                //Update scene and BVH
                Scene updatedScene {
                    .type = scene.type,
                    .meshes = meshes,
                    .spheres = spheres,
                    .lights = scene.lights

                };
                BVH bvh = BVH(updatedScene, features);

                RenderState updatedState = {
                    .scene = updatedScene,
                    .features = state.features,
                    .bvh = bvh,
                    .sampler = state.sampler
                };

                //Generate scene rays
                auto rays = generatePixelRays(updatedState, camera, { x, y }, screen.resolution());
                sceneLight += renderRays(updatedState, rays);
            }

            sceneLight /= (float) features.extra.numBlurSamples;
            screen.setPixel(x, y, sceneLight);
        }
    }

}



// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    // ...
}


// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...
}

// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    // DEPENDS ON TEXTURE MAPPING
    if (state.features.extra.enableEnvironmentMap && state.features.enableTextureMapping) {
        const std::shared_ptr<Image> env = state.scene.environment;
        const glm::vec3& n = glm::normalize(ray.direction);

        // Sources:
        // Marschner, S., & Shirley, P. (2015). Fundamentals of computer graphics (Fourth). CRC Press, Taylor & Francis Group. Retrieved October 30, 2023, from https://learning-oreilly-com.tudelft.idm.oclc.org/library/view/fundamentals-of-computer/9781482229417/K22616_C011.xhtml#:-:text=Spherical%20Coordinates,line%20of%20latitude.

        // Section 11.2.1 - Spherical Coordinates, from the book
        const float lambda = atan2(n.z, n.x);
        const float theta = glm::acos(n.y);

        // Map angles to texture coordinates
        const float x = (lambda + glm::pi<float>()) / glm::two_pi<float>();
        const float y = (glm::pi<float>() - theta) / glm::pi<float>();

        const glm::vec2 tex { x, y };

        // Use bilinear interpolation if available
        if (state.features.enableBilinearTextureFiltering) {
            return sampleTextureBilinear(*env, tex);
        }
        return sampleTextureNearest(*env, tex);
    }
    return glm::vec3(0.f);
}


// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    using Primitive = BVH::Primitive;

    return 0; // This is clearly not the solution
}