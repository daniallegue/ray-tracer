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


constexpr long double factorial(size_t n)
{
    if (n < 2)
        return 1;

    long double result = 1;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

constexpr float combination(size_t k, size_t n)
{
    //TODO use iterative nCr function
    return static_cast<float>(factorial(n) / (factorial(k) * factorial(n - k)));
}

std::vector<float> gaussianFilterVector(uint32_t n)
{
    std::vector<float> result(n);
    float sum = 0.0f;
    for (size_t i = 0; i < n; i++) {
        const float a = combination(i + 1, n);
        result[i] = a;
        sum += a;
    }
    for (size_t i = 0; i < n; i++) {
        result[i] = result[i] / sum;
    }

    return result;
}

std::vector<float> gaussian_filter (0);

glm::vec3 applyFilter(const std::vector<float>& filter, const uint32_t axis, const Screen& image, const std::vector<glm::vec3>& pixels, const int x, const int y)
{
    const size_t filtersize = filter.size();
    const size_t filterradius = filtersize / 2;

    if (axis > 1)
        return {};

    glm::vec3 total = {};

    if (axis == 0) {
        const int width = image.resolution().x;
        for (int i = 0; i < filtersize; i++) {
            const size_t index = std::clamp<size_t>(x + i - filterradius, 0, width - 1);
            total += pixels[image.indexAt(index, y)] * filter[i];
        }
    } else if (axis == 1) {
        const int height = image.resolution().y;
        for (int i = 0; i < filtersize; i++) {
            const size_t index = std::clamp<size_t>(y + i - filterradius, 0, height - 1);
            total += pixels[image.indexAt(x, index)] * filter[i];
        }
    }

    return total;
}

// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    if (bloom_filter_size != gaussian_filter.size()) {
        // Create new gaussian filter
        gaussian_filter = gaussianFilterVector(bloom_filter_size);
    }

    const std::vector<glm::vec3>& pixels = image.pixels();
    std::vector<glm::vec3> buffer0(image.pixels().size());
    std::vector<glm::vec3> buffer1(image.pixels().size());
    const int width = image.resolution().x;
    const int height = image.resolution().y;

    // Filter pixels by threshold and store in first buffer
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            const glm::vec3 v = pixels[image.indexAt(x, y)];
            if (v.x > bloom_threshold || v.y > bloom_threshold || v.z > bloom_threshold) {
                buffer0[image.indexAt(x, y)] = v;
            }
        }
    }

    // Apply horizontal blur filter and store in second buffer
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            buffer1[image.indexAt(x, y)] = applyFilter(gaussian_filter, 0, image, buffer0, x, y);
        }
    }

    // Apply vertical blur filter and add blurred result to the original image
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            glm::vec3 v = pixels[image.indexAt(x, y)] + applyFilter(gaussian_filter, 1, image, buffer1, x, y);
            image.setPixel(x, y, glm::clamp(v, 0.0f, 1.0f));
        }
    }
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