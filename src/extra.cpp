#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>

//own recursive method to calculate the resulting color
//it uses renderRay and then I calculate my own color
glm::vec3 rayTraceDepthOfField(RenderState& state, Ray ray, float apertureSize, float focalLength, int rayDepth)
{

    // calculate all other components
    auto L = renderRay(state, ray, rayDepth);

    // calculate the point that will be focused
    glm::vec3 intersection = ray.origin + ray.t * ray.direction;
    glm::vec3 focalPoint = intersection + focalLength * ray.direction;

    //randomly generated point within camera aperture
    glm::vec3 randomAperturePoint(0.0f);
    randomAperturePoint.x = (static_cast<float>(rand() % RAND_MAX) / RAND_MAX) * apertureSize - apertureSize / 2.0f;
    randomAperturePoint.y = (static_cast<float>(rand() % RAND_MAX) / RAND_MAX) * apertureSize - apertureSize / 2.0f;
    randomAperturePoint.z = (static_cast<float>(rand() % RAND_MAX) / RAND_MAX) * apertureSize - apertureSize / 2.0f;

    //this point will represent a blur point since shifted from origin
    glm::vec3 randomPointOnLens =ray.origin + randomAperturePoint;
    glm::vec3 directionOfRandomPoint = glm::normalize(intersection - randomAperturePoint);

    
    //this ray simulates where the camera will be focused, essentially from randomAperturePoint in the direction of apertureToFocal
    Ray r(randomPointOnLens, directionOfRandomPoint);

    rayDepth += 1;
    L += rayTraceDepthOfField(state, r, apertureSize, focalLength, rayDepth);
    return L;
}

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
    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            RenderState state = {
                .scene = scene,
                .features = features,
                .bvh = bvh,
                .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
            };
            auto rays = generatePixelRays(state, camera, { x, y }, screen.resolution());

            glm::vec3 diffuseColor { 0.f };
            for (const auto& ray : rays) {
                diffuseColor += rayTraceDepthOfField(state,ray, static_cast<float>(features.extra.apertureSize), static_cast<float>(features.extra.focalLength),0);
            }
            diffuseColor /= static_cast<float>(rays.size());
            screen.setPixel(x, y, diffuseColor);
        }
    }
    // ...
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
    auto numSamples = state.features.extra.numGlossySamples;
    if (hitInfo.material.ks == glm::vec3(0.0f)) {
        return;
    }
    glm::vec3 intersection = ray.origin + ray.t * ray.direction;
    glm::vec3 reflection = glm::normalize(glm::reflect(ray.direction, hitInfo.normal));

    glm::vec3 glossyColor(0.0f, 0.0f, 0.0f);
    std::vector<Ray> glossyRays;

    //generate numSamples Glossy Rays
    for (int x = 0; x < numSamples; x++) {
        //Generates and angle between [0,2pi]
        float angle = static_cast<float>(2.0 * glm::pi<float>() * rand() / static_cast<float>(RAND_MAX));

        //Generate a radius between [0,1]
        float radius = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);

        //taking the corresponding disk amplitude
        float a = hitInfo.material.shininess / 256.0f;
        radius = a * sqrt(radius);
        
        //Convert polar coordinates to Cartesian coordinates
        float r1 = radius * std::cos(angle);
        float r2 = radius * std::sin(angle);

        //disturb reflection
        glm::vec3 u = hitInfo.normal;
        glm::vec3 v = glm::cross(u, reflection);
        glm::vec3 pertubedReflection = glm::normalize(reflection + u*r1 + v*r2);

        Ray glossyRay(intersection + pertubedReflection * 10.0f * std::numeric_limits<float>::epsilon(), pertubedReflection);
        
        glossyRays.push_back(glossyRay);
    }

    //calculate the final Color for each glossy Ray and taking its corresponding weight ->  /(x+1)
    for (float x = 0.0f; x < glossyRays.size(); x++) {
        hitColor = hitColor * x / (x + 1) + (renderRay(state, glossyRays[x], rayDepth + 1) * hitInfo.material.ks) / (x + 1);
    }
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        return glm::vec3(0.f);
    } else {
        return glm::vec3(0.f);
    }
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