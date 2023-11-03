#include "extra.h"
#include "bvh.h"
#include "intersect.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include "texture.h"
#include "draw.h"
#include <framework/trackball.h>
#include <numeric>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!

//Citations: https://pathtracing.home.blog/depth-of-field/

void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }
    float apertureSize = static_cast<float>(features.extra.apertureSize);
    float focalLength = static_cast<float>(features.extra.focalLength);
    float numSamples = static_cast<float>(features.extra.numDOFSamples);
    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            RenderState state = {
                .scene = scene,
                .features = features,
                .bvh = bvh,
                .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
            };
            auto rays = generatePixelRays(state, camera, { x, y }, screen.resolution());
            std::vector<Ray> newRays;  
            glm::vec3 newColor { 0.f };
            for (float samples = 0.0f; samples < numSamples; samples++) {
                for (const auto& ray : rays) {

                    // randomly generated point within aperture
                    glm::vec3 randomAperturePoint(0.0f);
                    randomAperturePoint.x = (static_cast<float>(rand() % RAND_MAX) / RAND_MAX) * apertureSize - apertureSize / 2.0f;
                    randomAperturePoint.y = (static_cast<float>(rand() % RAND_MAX) / RAND_MAX) * apertureSize - apertureSize / 2.0f;
                    randomAperturePoint.z = (static_cast<float>(rand() % RAND_MAX) / RAND_MAX) * apertureSize - apertureSize / 2.0f;

                    // New intersection refering to focal length
                    glm::vec3 newIntersection = ray.origin + focalLength * ray.direction;

                    // new origin with random aperture
                    glm::vec3 newOrigin = ray.origin + randomAperturePoint;

                    // new direction from newOrigin to newIntersection
                    glm::vec3 newDirection = glm::normalize(newIntersection - newOrigin);

                    // find value of T
                    float newT = (newIntersection.x - newOrigin.x) / newDirection.x;

                    // new Ray
                    Ray newRay = { newOrigin, newDirection };

                    newRays.push_back(newRay);
                }
                // render all set of new rays
                newColor += renderRays(state, newRays);
            }
            
            screen.setPixel(x, y, newColor/numSamples);
        }
    }
    // ...
    //diffuseColor += rayTraceDepthOfField(state, ray,, , 0);
}

//Define function for spline transformation with a cubic B�zier curve
glm::mat4 cubicBezierTransformation(const Features& features, glm::vec3 pos, float time)
{
    float u = 1.0f - time;
    float u2 = glm::pow(u, 2);
    float u3 = glm::pow(u, 3);
    float t2 = glm::pow(time, 2);
    float t3 = glm::pow(time, 3);


    glm::vec3 p0 = (glm::vec3(0, 0, 0) * 1.02f) + pos;
    glm::vec3 p1 = (glm::vec3(1, 2, 2) * 1.02f) + pos;
    glm::vec3 p2 = (glm::vec3(1, 2, 2) * 1.02f) + pos;
    glm::vec3 p3 = (glm::vec3(3, 1, 0) * 1.02f) + pos;
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
    //Resources:
    //1. Fundamentals of Computer Graphics, 4th Edition, Marschner & Shirley, Section 13.4.5
    //2. https://en.wikipedia.org/wiki/Motion_blur
    //3. https://en.wikipedia.org/wiki/B%C3%A9zier_curve
    //4. https://math.stackexchange.com/questions/867153/what-is-the-parametric-function-of-the-new-bezier-curve
    //5. Answers EWI

    if (!features.extra.enableMotionBlur) {
        return;
    }
#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
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
                //Generate sample for iteration in the range [0, 1]
                float time = state.sampler.next_1d();

                //Apply to meshes
                std::vector<Mesh> meshes;
                for (int i = 0; i < scene.meshes.size(); i++) {
                    Mesh mesh = scene.meshes[i];
                    std::vector<Vertex> updatedVertices;
                    std::vector<Vertex> vertices = mesh.vertices;
                    //Apply transformation to each vertex
                    for (int u = 0; u < vertices.size(); u++) {
                        Vertex v = vertices[u];
                        glm::vec3 pos = v.position;
                        glm::mat4 transform = cubicBezierTransformation(features, pos, time);
                        glm::vec4 newPos = transform * glm::vec4 { pos.x, pos.y, pos.z, 1 };
                        //3-Dimensional vector
                        glm::vec3 posScaled = { newPos[0] / newPos[3], newPos[1] / newPos[3], newPos[2] / newPos[3] };
                        // DrawLine for visuall debugging purposes
                        drawSegment(pos, posScaled, glm::vec3 { 1.0f, 1.0f, 0.0f });

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
                    //Apply transformation to the sphere's center
                    Sphere sphere = scene.spheres[i];
                    glm::vec3 pos = sphere.center;
                    glm::mat4 transform = cubicBezierTransformation(features, pos, time);
                    glm::vec4 newPos = transform * glm::vec4 { pos.x, pos.y, pos.z, 1 };
                    // 3-Dimensional vector
                    glm::vec3 posScaled = { newPos[0] / newPos[3], newPos[1] / newPos[3], newPos[2] / newPos[3] };
                    // DrawLine for visuall debugging purposes
                    glm::vec3 direction = posScaled - pos;
                    drawLine(pos, direction, glm::vec3 { 0.0f, 1.0f, 0.0f });
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

            //Average values
            sceneLight /= (float) features.extra.numBlurSamples;
            screen.setPixel(x, y, sceneLight);
        }
    }

}



// TODO; Extra feature

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
    // TODO use iterative nCr function
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

std::vector<float> gaussian_filter(0);

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

//Citations: Book Fundamental of Computer Graphics 4th edition by Peter Shirley Pages 333,334

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

        //taking the corresponding disk amplitude, in this case I chose:
        float a = hitInfo.material.shininess / 64.0f;
        radius = a * sqrt(radius);
        
        //Convert polar coordinates to Cartesian coordinates
        float r1 = radius * std::cos(angle);
        float r2 = radius * std::sin(angle);

        //created an orthonormal basis
        glm::vec3 u = hitInfo.normal;
        glm::vec3 v = glm::cross(u, reflection);

        //sampled on that basis based on the disk
        glm::vec3 pertubedReflection = glm::normalize(reflection + u*r1 + v*r2);

        Ray glossyRay(intersection + pertubedReflection * 10.0f * std::numeric_limits<float>::epsilon(), pertubedReflection);
        
        glossyRays.push_back(glossyRay);
    }

    //calculate the final Color for each glossy Ray and taking its corresponding weight ->  /(x+1)
    for (float x = 0.0f; x < glossyRays.size(); x++) {
        hitColor = hitColor * x / (x + 1) + (renderRay(state, glossyRays[x], rayDepth + 1) * hitInfo.material.ks) / (x + 1);
    }
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
        if (!env) {
            // Environment is missing (should not happen)
            assert(env);
            return {};
        }

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

#pragma region SAH+Binning

/// <summary>
///  Calculates the surface area of an AABB
/// </summary>
/// <param name="aabb"></param>
/// <returns>Surface area of aabb</returns>
float computeSurfaceAreaAABB(const AxisAlignedBox& aabb)
{
    const glm::vec3 v = aabb.upper - aabb.lower;

    // IMPORTANT! The data structures lecture slides mention using surface area. 
    // This seems weird as this is in 3D. Even the 2019-2020 lecture recording about acceleration data structures (1:26:45) mentions the 3D surface area is the volume.
    // However, this does not add up for a SAH for a BVH for ray tracing, as a bounding box that is completely squashed with 0 volume, can still be perfectly intersected by a ray.
    // Therefore I have chosen to use the actual surface area of the 6 sides of the bounding box.

    return 2 * (v.x * v.y + v.y * v.z + v.z * v.x);
}

AxisAlignedBox computeSpanAABB(const std::span<const AxisAlignedBox> aabbs)
{
    if (aabbs.size() <= 0)
        return AxisAlignedBox {};

    AxisAlignedBox acc = aabbs[0];
    for (const AxisAlignedBox& aabb : aabbs)
        acc = UnionAABB(acc, aabb);
    return acc;
}

float splitcost(const std::span<AxisAlignedBox>& aabbs, uint32_t splitindex, const AxisAlignedBox& aabb)
{
    const size_t num_a = splitindex;
    const size_t num_b = aabbs.size() - splitindex;

    if (num_a == 0 || num_b == 0) {
        // If this occurs, we need to make sure it doesn't get split,
        // otherwise the child node also doesn't get split and will
        // continue on forever, resulting in stackoverflow.
        return std::numeric_limits<float>::infinity();
    }

    // CITE: Data structures lecture slides page 144
    //   aabb0 corresponds to A
    //   aabb1 corresponds to B
    const AxisAlignedBox aabb0 = computeSpanAABB(aabbs.subspan(0, num_a));
    const AxisAlignedBox aabb1 = computeSpanAABB(aabbs.subspan(num_a, num_b));

    // CITE: Data structures lecture slides page 148
    // Calculate naive probability of needing to intersect children using surface area (volume).
    const float surfacearea = computeSurfaceAreaAABB(aabb);
    const float p_a = computeSurfaceAreaAABB(aabb0) / surfacearea;
    const float p_b = computeSurfaceAreaAABB(aabb1) / surfacearea;

    // Must always be true. (Only executed in debug mode)
    assert(p_a <= 1 && p_a >= 0);
    assert(p_b <= 1 && p_b >= 0);

    // CITE: Data structures lecture slides page 146
    // Using constant costs for traversal and intersection (0 and 1 respectively)
    const float cost = p_a * num_a + p_b * num_b;
    return cost;
}

glm::vec3 computeSplitPoint(const AxisAlignedBox& aabb, const uint32_t axis, const uint32_t num_bins, const uint32_t bin)
{
    glm::vec3 v {};
    v[axis] = (aabb.upper - aabb.lower)[axis];
    return aabb.lower + v * (static_cast<float>(bin + 1) / num_bins);
}

int computeSplitIndex(const std::span<float>& centroids, const uint32_t axis, const glm::vec3& splitpoint)
{
    const float boundary = splitpoint[axis];
    for (int i = 0; i < centroids.size(); i++) {
        // Make sure list is sorted (Only executed in debug mode)
        float centroidf = centroids[i];
        if (i < centroids.size() - 1) {
            assert(centroidf <= centroids[i + 1]);
        }

        if (centroidf > boundary) {
            return i;
        }
    }
    // Fallback, split in the middle
    return centroids.size() / 2;
}

std::vector<glm::vec3> computeSplitPoints(const AxisAlignedBox& aabb, const uint32_t axis, const uint32_t num_bins)
{
    std::vector<glm::vec3> splitpoints(num_bins - 1);
    for (int i = 0; i < num_bins - 1; i++) {
        const glm::vec3 point = computeSplitPoint(aabb, axis, num_bins, i);
        splitpoints[i] = point;
    }
    return splitpoints;
}

std::vector<uint32_t> computeSplitIndexes(const std::span<float>& axialcentroids, const std::vector<glm::vec3>& splitpoints, const uint32_t axis)
{
    std::vector<uint32_t> splitindexes(splitpoints.size());
    std::transform(splitpoints.begin(), splitpoints.end(), splitindexes.begin(), [axialcentroids, axis](const glm::vec3& splitpoint) {
        return computeSplitIndex(axialcentroids, axis, splitpoint);
    });
    return splitindexes;
}

std::vector<BVH::SplitPlane> calculateSplitPlanes(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    if (axis > 2)
        return {};

    std::ranges::sort(primitives, [axis](const BVHInterface::Primitive& a, const BVHInterface::Primitive& b) {
        return computePrimitiveCentroid(a)[axis] < computePrimitiveCentroid(b)[axis];
    });

    // CITE: Data structures lecture slides page 152
    // Calculate number of bins based on number of primitives with a logarithmic function such that "#bins << #primitives"
    const uint32_t numbins = num_bins(primitives.size());

    // Compute axis components of centroids for all primitives
    std::vector<float> centroids(primitives.size());
    std::ranges::transform(primitives, centroids.begin(), [axis](BVH::Primitive p) {
        return computePrimitiveCentroid(p)[axis];
    });

    // Compute aabbs for all primitives
    std::vector<AxisAlignedBox> aabbs(primitives.size());
    std::ranges::transform(primitives, aabbs.begin(), [](BVH::Primitive p) {
        return computePrimitiveAABB(p);
    });
    std::span<AxisAlignedBox> aabbspan { aabbs };

    const std::vector<glm::vec3> splitpoints = computeSplitPoints(aabb, axis, numbins);
    const std::vector<uint32_t> splitindexes = computeSplitIndexes(centroids, splitpoints, axis);

    // CITE: Data structures lecture slides page 152
    // Calculate the split costs for all split indexes
    std::vector<float> costs(splitindexes.size());
    std::transform(splitindexes.begin(), splitindexes.end(), costs.begin(), [aabbspan, aabb](uint32_t splitindex) {
        return splitcost(aabbspan, splitindex, aabb);
    });

    // Create split planes
    std::vector<BVH::SplitPlane> splitplanes(splitindexes.size());
    for (int i = 0; i < splitindexes.size(); i++) {
        const uint32_t splitindex = splitindexes[i];
        const glm::vec3 splitpoint = splitpoints[i];
        const float cost = costs[i];
        splitplanes[i] = BVH::SplitPlane { splitpoint, axis, splitindex, cost, glm::vec4 {} };
    }
    return splitplanes;
}

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
    using SplitPlane = BVH::SplitPlane;
    // CITE: Data structures lecture slides page 152
    // Calculate "predefined splitplanes" and "compute cost for each plane"
    std::vector<SplitPlane> splitplanes = calculateSplitPlanes(aabb, axis, primitives);

    // Remove any split indexes that don't create an actual split to prevent infinite recursion
    const size_t len = primitives.size();
    std::ranges::remove_if(splitplanes, [len](const BVH::SplitPlane &a){ 
        return a.index <= 0 || a.index >= len - 1;
    });

    if (splitplanes.empty()) {
        // Fallback
        return splitPrimitivesByMedian(aabb, axis, primitives);
    }

    // Find splitplane with minimal cost (pointer will never be nullptr, because list is not empty)
    const SplitPlane splitplane = *std::ranges::min_element(splitplanes, [](const SplitPlane &a, const SplitPlane &b) {
        return a.cost < b.cost;
    });
    return splitplane.index;
}

#pragma endregion
