#pragma once

#include "fwd.h"
#include "bvh.h"
#include "bvh_interface.h"
#include "render.h"
#include "scene.h"
#include "screen.h"

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of depth of field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen);

//Provides a spline matrix transformation given the time and intial position
glm::mat4 cubicBezierTransformation(const Features& features, glm::vec3 pos, float time);
// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// allowing objects to move during a render, and visualize the appearance of movement.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen);

inline int bloom_filter_size = 41;
inline float bloom_threshold = 0.9f;
inline float bloom_scalar = 1.5f;

// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& screen);

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// For a description of the method's arguments, refer to 'extra.cpp'
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth);

// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You may have to add support for environment maps
// to the Scene object, and provide a scene with the right data for this.
// For a description of the method's arguments, refer to 'extra.cpp'
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray);

/// <summary>
///  Calculates a appropriate number of bins according to the number of nodes.
/// </summary>
/// <param name="numnodes"></param>
/// <returns>Number of bins</returns>
constexpr uint32_t num_bins(uint32_t numnodes)
{
    return std::max(2u, static_cast<uint32_t>(std::ceil(std::log2(numnodes))));
}

/// <summary>
///  Calculate the cost to split the given set of nodes at the splitindex
/// </summary>
/// <param name="aabbs">Span of <see cref="::AxisAlignedBox"/></param>
/// <param name="splitindex"></param>
/// <param name="aabb"></param>
/// <returns></returns>
float splitcost(const std::span<AxisAlignedBox>& aabbs, uint32_t splitindex, const AxisAlignedBox& aabb);

/// <summary>
///  Computes num_bins - 1 split points that lie on the split planes
/// </summary>
/// <param name="aabb">Bounding box to generate planes in</param>
/// <param name="axis"></param>
/// <param name="num_bins"></param>
/// <returns>Points on the split planes inside the bounding box</returns>
std::vector<glm::vec3> computeSplitPoints(const AxisAlignedBox& aabb, const uint32_t axis, const uint32_t num_bins);

/// <summary>
///  Converts the split points into a list of split indexes.
/// </summary>
/// <param name="primitives">List of the axis components of the primitives to base the index upon.</param>
/// <param name="splitpoints"></param>
/// <param name="axis"></param>
/// <returns></returns>
std::vector<uint32_t> computeSplitIndexes(const std::span<float>& axialcentroids, const std::vector<glm::vec3>& splitpoints, const uint32_t axis);

/// <summary>
///  Calculates split planes based on a list of primitives
/// </summary>
/// <param name="aabb">Bounding box around all primitives</param>
/// <param name="axis"></param>
/// <param name="primitives"></param>
/// <returns>Splitplanes</returns>
std::vector<BVH::SplitPlane> calculateSplitPlanes(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives);

// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// For a description of the method's arguments, refer to 'bounding_volume_hierarchy.cpp'
// NOTE: this method is unit-tested, so do not change the function signature.
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives);