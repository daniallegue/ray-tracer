#include "bvh.h"
#include "draw.h"
#include "extra.h"
#include "interpolate.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "texture.h"
#include <algorithm>
#include <bit>
#include <chrono>
#include <framework/opengl_includes.h>
#include <iostream>
#include <queue>
#include <stack>

// Helper method to fill in hitInfo object. This can be safely ignored (or extended).
// Note: many of the functions in this helper tie in to standard/extra features you will have
// to implement separately, see interpolate.h/.cpp for these parts of the project
void updateHitInfo(RenderState& state, const BVHInterface::Primitive& primitive, const Ray& ray, HitInfo& hitInfo)
{
    const auto& [v0, v1, v2] = std::tie(primitive.v0, primitive.v1, primitive.v2);
    const auto& mesh = state.scene.meshes[primitive.meshID];
    const auto n = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
    const auto p = ray.origin + ray.t * ray.direction;

    // First, fill in default data, unrelated to separate features
    hitInfo.material = mesh.material;
    hitInfo.normal = n;
    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, p);

    // Next, if `features.enableNormalMapping` is true, generate smoothly interpolated vertex normals
    if (state.features.enableNormalInterp) {
        drawLine(v0.position, v0.normal, glm::vec3(0.5f, 0.5f, 0.5f));
        drawLine(v1.position, v1.normal, glm::vec3(0.5f, 0.5f, 0.5f));
        drawLine(v2.position, v2.normal, glm::vec3(0.5f, 0.5f, 0.5f));
        glm::vec3 normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
        drawLine(p, normal, glm::vec3(0.5f, 0.5f, 0.5f));
        hitInfo.normal = normal;
    }

    // Next, if `features.enableTextureMapping` is true, generate smoothly interpolated vertex uvs
    if (state.features.enableTextureMapping) {
        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
    }

    // Finally, catch flipped normals
    if (glm::dot(ray.direction, n) > 0.0f) {
        hitInfo.normal = -hitInfo.normal;
    }
}

// BVH constructor; can be safely ignored. You should not have to touch this
// NOTE: this constructor is tested, so do not change the function signature.
BVH::BVH(const Scene& scene, const Features& features)
{
#ifndef NDEBUG
    // Store start of bvh build for timing
    using clock = std::chrono::high_resolution_clock;
    const auto start = clock::now();
#endif

    // Count the total nr. of triangles in the scene
    size_t numTriangles = 0;
    for (const auto& mesh : scene.meshes)
        numTriangles += mesh.triangles.size();

    // Given the input scene, gather all triangles over which to build the BVH as a list of Primitives
    std::vector<Primitive> primitives;
    primitives.reserve(numTriangles);
    for (uint32_t meshID = 0; meshID < scene.meshes.size(); meshID++) {
        const auto& mesh = scene.meshes[meshID];
        for (const auto& triangle : mesh.triangles) {
            primitives.push_back(Primitive {
                .meshID = meshID,
                .v0 = mesh.vertices[triangle.x],
                .v1 = mesh.vertices[triangle.y],
                .v2 = mesh.vertices[triangle.z] });
        }
    }

    // Tell underlying vectors how large they should approximately be
    m_primitives.reserve(numTriangles);
    m_nodes.reserve(numTriangles + 1);

    // Recursively build BVH structure; this is where your implementation comes in
    m_nodes.emplace_back(); // Create root node
    m_nodes.emplace_back(); // Create dummy node s.t. children are allocated on the same cache line
    if (debugBuildBVH) {
        // Build initial giant leaf
        m_nodes[RootIndex] = buildLeafData(scene, features, computeSpanAABB(primitives), primitives);
    } else {
        // Build normal
        buildRecursive(scene, features, primitives, RootIndex);
    }

    // Fill in boilerplate data
    buildNumLevels();
    buildNumLeaves();

#ifndef NDEBUG
    // Output end of bvh build for timing
    const auto end = clock::now();
    std::cout << "BVH construction time: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms" << std::endl;
#endif
}

// BVH helper method; allocates a new node and returns its index
// You should not have to touch this
uint32_t BVH::nextNodeIdx()
{
    const auto idx = static_cast<uint32_t>(m_nodes.size());
    m_nodes.emplace_back();
    return idx;
}

// Given a BVH triangle, compute an axis-aligned bounding box around the primitive
// - primitive; a single triangle to be stored in the BVH
// - return;    an axis-aligned bounding box around the triangle
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computePrimitiveAABB(const BVHInterface::Primitive primitive)
{
    const std::array vertices { primitive.v1.position, primitive.v2.position };
    glm::vec3 lower = primitive.v0.position;
    glm::vec3 upper = primitive.v0.position;
    for (const glm::vec3& vertex : vertices) {
        if (vertex.x < lower.x)
            lower.x = vertex.x;
        if (vertex.x > upper.x)
            upper.x = vertex.x;

        if (vertex.y < lower.y)
            lower.y = vertex.y;
        if (vertex.y > upper.y)
            upper.y = vertex.y;

        if (vertex.z < lower.z)
            lower.z = vertex.z;
        if (vertex.z > upper.z)
            upper.z = vertex.z;
    }

    return { lower, upper };
}

// Given a range of BVH triangles, compute an axis-aligned bounding box around the range.
// - primitive; a contiguous range of triangles to be stored in the BVH
// - return;    a single axis-aligned bounding box around the entire set of triangles
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computeSpanAABB(std::span<const BVHInterface::Primitive> primitives)
{
    const size_t size = primitives.size();
    if (size == 0)
        return {};
    AxisAlignedBox aabb = computePrimitiveAABB(primitives[0]);
    for (size_t i = 1; i < size; i++) {
        aabb = UnionAABB(aabb, computePrimitiveAABB(primitives[i]));
    }
    return aabb;
}

AxisAlignedBox UnionAABB(const AxisAlignedBox& a, const AxisAlignedBox& b)
{
    const std::array vertices { a.upper, b.lower, b.upper };
    glm::vec3 lower = a.lower;
    glm::vec3 upper = a.lower;
    for (const glm::vec3& vertex : vertices) {
        if (vertex.x < lower.x)
            lower.x = vertex.x;
        if (vertex.x > upper.x)
            upper.x = vertex.x;

        if (vertex.y < lower.y)
            lower.y = vertex.y;
        if (vertex.y > upper.y)
            upper.y = vertex.y;

        if (vertex.z < lower.z)
            lower.z = vertex.z;
        if (vertex.z > upper.z)
            upper.z = vertex.z;
    }

    return { lower, upper };
}

// Given a BVH triangle, compute the geometric centroid of the triangle
// - primitive; a single triangle to be stored in the BVH
// - return;    the geometric centroid of the triangle's vertices
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePrimitiveCentroid(const BVHInterface::Primitive primitive)
{
    return (primitive.v0.position + primitive.v1.position + primitive.v2.position) * glm::third<float>();
}

// Given an axis-aligned bounding box, compute the longest axis; x = 0, y = 1, z = 2.
// - aabb;   the input axis-aligned bounding box
// - return; 0 for the x-axis, 1 for the y-axis, 2 for the z-axis
//           if several axes are equal in length, simply return the first of these
// This method is unit-tested, so do not change the function signature.
uint32_t computeAABBLongestAxis(const AxisAlignedBox& aabb)
{
    const glm::vec3 v = aabb.upper - aabb.lower;
    if (v.x >= v.y && v.x >= v.z)
        return 0;
    if (v.y >= v.z)
        return 1;
    return 2;
}

// Given a range of BVH triangles, sort these along a specified axis based on their geometric centroid.
// Then, find and return the split index in the range, such that the subrange containing the first element
// of the list is at least as big as the other, and both differ at most by one element in size.
// Hint: you should probably reuse `computePrimitiveCentroid()`
// - aabb;       the axis-aligned bounding box around the given triangle range
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires sorting/splitting along an axis
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesByMedian(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives)
{
    const size_t middle = (primitives.size() + 1) / 2;

    if (axis > 2)
        return middle;

    std::ranges::sort(primitives, [axis](const BVHInterface::Primitive& a, const BVHInterface::Primitive& b) {
        return computePrimitiveCentroid(a)[axis] < computePrimitiveCentroid(b)[axis];
    });

    return middle;
}

bool intersectRayWithAABB(const AxisAlignedBox& aabb, Ray &ray)
{
    // CITE: Ray Tracing Lecture Slides page 81 
    float oneoverx = 1.0f / ray.direction.x;
    float oneovery = 1.0f / ray.direction.y;
    float oneoverz = 1.0f / ray.direction.z;

    float txmin = (aabb.lower.x - ray.origin.x) * oneoverx;
    float txmax = (aabb.upper.x - ray.origin.x) * oneoverx;
    float tymin = (aabb.lower.y - ray.origin.y) * oneovery;
    float tymax = (aabb.upper.y - ray.origin.y) * oneovery;
    float tzmin = (aabb.lower.z - ray.origin.z) * oneoverz;
    float tzmax = (aabb.upper.z - ray.origin.z) * oneoverz;

    float tin = std::max(std::min(txmin, txmax), std::max(std::min(tymin, tymax), std::min(tzmin, tzmax)));
    float tout = std::min(std::max(txmin, txmax), std::min(std::max(tymin, tymax), std::max(tzmin, tzmax)));

    if (tin > tout || tout < 0)
        return false;

    return true;
}

// Hierarchy traversal routine; called by the BVH's intersect(),
// you must implement this method and implement it carefully!
//
// If `features.enableAccelStructure` is not enabled, the method should just iterate the BVH's
// underlying primitives (or the scene's geometry). The default imlpementation already does this.
// You will have to implement the part which actually traverses the BVH for a faster intersect,
// given that `features.enableAccelStructure` is enabled.
//
// This method returns `true` if geometry was hit, and `false` otherwise. On first/closest hit, the
// distance `t` in the `ray` object is updated, and information is updated in the `hitInfo` object.
//
// - state;    the active scene, and a user-specified feature config object, encapsulated
// - bvh;      the actual bvh which should be traversed for faster intersection
// - ray;      the ray intersecting the scene's geometry
// - hitInfo;  the return object, with info regarding the hit geometry
// - return;   boolean, if geometry was hit or not
//
// This method is unit-tested, so do not change the function signature.
bool intersectRayWithBVH(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo)
{
    // Relevant data in the constructed BVH
    const std::span<const BVHInterface::Node> nodes = bvh.nodes();
    const std::span<const BVHInterface::Primitive> primitives = bvh.primitives();

    // Return value
    bool is_hit = false;

    if (state.features.enableAccelStructure) {
        //
        // Some hints (refer to bvh_interface.h either way). BVH nodes are packed, so the
        // data is not easily extracted. Helper methods are available, however:
        // - For a given node, you can test if the node is a leaf with `node.isLeaf()`.
        // - If the node is not a leaf, you can obtain the left/right children with `node.leftChild()` etc.
        // - If the node is a leaf, you can obtain the offset to and nr. of primitives in the bvh's list
        //   of underlying primitives with `node.primitiveOffset()` and `node.primitiveCount()`
        //
        // In short, you will have to step down the bvh, node by node, and intersect your ray
        // with the node's AABB. If this intersection passes, you should:
        // - if the node is a leaf, intersect with the leaf's primitives
        // - if the node is not a leaf, test the left and right children as well!
        //
        // Note that it is entirely possible for a ray to hit a leaf node, but not its primitives,
        // and it is likewise possible for a ray to hit both children of a node.

        // Using stack-based Depth-first traversal
        if (!nodes.empty()) {
            std::stack<uint32_t> nodestack {};
            nodestack.push(BVH::RootIndex);

            while (!nodestack.empty()) {
                const BVHInterface::Node& node = nodes[nodestack.top()];
                nodestack.pop();

                if (intersectRayWithAABB(node.aabb, ray)) {
                    if (node.isLeaf()) {

                        const uint32_t count = node.primitiveCount();
                        const uint32_t offset = node.primitiveOffset();

                        for (uint32_t i = 0; i < count; i++) {
                            const BVHInterface::Primitive& primitive = primitives[offset + i];
                            const auto& [v0, v1, v2] = std::tie(primitive.v0, primitive.v1, primitive.v2);
                            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                                updateHitInfo(state, primitive, ray, hitInfo);
                                is_hit = true;
                            }
                        }
                    } else {
                        nodestack.push(node.leftChild());
                        nodestack.push(node.rightChild());
                    }
                } 
            }
        }
    } else {
        // Naive implementation; simply iterates over all primitives
        for (const auto& prim : primitives) {
            const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                updateHitInfo(state, prim, ray, hitInfo);
                is_hit = true;
            }
        }
    }

    // Intersect with spheres.
    for (const auto& sphere : state.scene.spheres)
        is_hit |= intersectRayWithShape(sphere, ray, hitInfo);

    return is_hit;
}

// Leaf construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and a range of triangles, generate a valid leaf object
// and store the triangles in the `m_primitives` vector.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;      the active scene
// - features;   the user-specified features object
// - aabb;       the axis-aligned bounding box around the primitives beneath this leaf
// - primitives; the range of triangles to be stored for this leaf
BVH::Node BVH::buildLeafData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, std::span<Primitive> primitives)
{
    Node node;

    // When using the iterative BVH building, the m_primitives vector is already fully populated, so no need for copying the data into m_primitives since, they are already in the list
    // We only need to set the primitiveoffset
    if (primitives.begin()._Myptr >= m_primitives.begin()._Ptr && primitives.begin()._Myptr <= m_primitives.end()._Ptr) {
        node.data[0] = static_cast<uint32_t>(primitives.begin()._Myptr - m_primitives.begin()._Ptr) | Node::LeafBit;
    } else {
        node.data[0] = static_cast<uint32_t>(m_primitives.size()) | Node::LeafBit;
        // Copy the current set of primitives to the back of the primitives vector
        std::ranges::copy(primitives, std::back_inserter(m_primitives));
    }
    node.data[1] = static_cast<uint32_t>(primitives.size());
    node.aabb = aabb;

    return node;
}

// Node construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and left/right child indices, generate a valid node object.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;           the active scene
// - features;        the user-specified features object
// - aabb;            the axis-aligned bounding box around the primitives beneath this node
// - leftChildIndex;  the index of the node's left child in `m_nodes`
// - rightChildIndex; the index of the node's right child in `m_nodes`
BVH::Node BVH::buildNodeData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, uint32_t leftChildIndex, uint32_t rightChildIndex)
{
    Node node;

    node.data[0] = leftChildIndex & (~Node::LeafBit);
    node.data[1] = rightChildIndex & (~Node::LeafBit);
    node.aabb = aabb;

    return node;
}

// Hierarchy construction routine; called by the BVH's constructor,
// you must implement this method and implement it carefully!
//
// You should implement the other BVH standard features first, and this feature last, as you can reuse
// most of the other methods to assemble this part. There are detailed instructions inside the
// method which we recommend you follow.
//
// Arguments:
// - scene;      the active scene
// - features;   the user-specified features object
// - primitives; a range of triangles to be stored in the BVH
// - nodeIndex;  index of the node you are currently working on, this is already allocated
//
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildRecursive(const Scene& scene, const Features& features, std::span<Primitive> primitives, uint32_t nodeIndex)
{
    // WARNING: always use nodeIndex to index into the m_nodes array. never hold a reference/pointer,
    // because a push/emplace (in ANY recursive calls) might grow vectors, invalidating the pointers.

    // Compute the AABB of the current node.
    const AxisAlignedBox aabb = computeSpanAABB(primitives);

    // As a starting point, we provide an implementation which creates a single leaf, and stores
    // all triangles inside it. You should remove or comment this, and work on your own recursive
    // construction algorithm that implements the following steps. Make sure to reuse the methods
    // you have previously implemented to simplify this process.
    //
    // 1. Determine if the node should be a leaf, when the nr. of triangles is less or equal to 4
    //    (hint; use the `LeafSize` constant)
    // 2. If it is a leaf, fill in the leaf's data, and store its range of triangles in `m_primitives`
    // 3. If it is a node:
    //    3a. Split the range of triangles along the longest axis into left and right subspans,
    //        using either median or SAH-Binning based on the `Features` object
    //    3b. Allocate left/right child nodes
    //        (hint: use `nextNodeIdx()`)
    //    3c. Fill in the current node's data; aabb, left/right child indices
    //    3d. Recursively build left/right child nodes over their respective triangles
    //        (hint; use `std::span::subspan()` to split into left/right ranges)

    if (features.enableAccelStructure) {
        if (primitives.size() > BVH::LeafSize) {
            // Current node is an internal node

            size_t splitindex;
            const uint32_t axis = computeAABBLongestAxis(aabb);
            if (features.extra.enableBvhSahBinning) {
               
                std::vector<AxisAlignedBox> aabbs(primitives.size());

                // Consider all 3 axis
                float mincost = FLT_MAX;
                int32_t minaxis = -1;
                for (int i = 0; i < 3; i++) {
                    size_t idx = splitPrimitivesBySAHBin(aabb, i, primitives);
                    std::ranges::transform(primitives, aabbs.begin(), [](BVH::Primitive p) {
                        return computePrimitiveAABB(p);
                    });
                    float cost = splitcost(aabbs, idx, aabb);
                    if (cost < mincost) {
                        splitindex = idx;
                        mincost = cost;
                        minaxis = i;
                    }
                }
                if (minaxis == -1) {
                    // Fallback using median split
                    splitindex = splitPrimitivesByMedian(aabb, axis, primitives);
                } else {
                    // Sort using the best axis
                    std::ranges::sort(primitives, [minaxis](const BVHInterface::Primitive& a, const BVHInterface::Primitive& b) {
                        return computePrimitiveCentroid(a)[minaxis] < computePrimitiveCentroid(b)[minaxis];
                    });
                }
            } else {
                splitindex = splitPrimitivesByMedian(aabb, axis, primitives);
            }

            const auto leftprim = primitives.subspan(0, splitindex);
            const auto rightprim = primitives.subspan(splitindex, primitives.size() - splitindex);

            const uint32_t leftidx = nextNodeIdx();
            const uint32_t rightidx = nextNodeIdx();

            if (debugBuildBVH) {
                // Build iteratively
                const AxisAlignedBox leftaabb = computeSpanAABB(leftprim);
                const AxisAlignedBox rightaabb = computeSpanAABB(rightprim);
                m_nodes[leftidx] = buildLeafData(scene, features, leftaabb, leftprim);
                m_nodes[rightidx] = buildLeafData(scene, features, rightaabb, rightprim);
            } else {
                buildRecursive(scene, features, leftprim, leftidx);
                buildRecursive(scene, features, rightprim, rightidx);
            }

            m_nodes[nodeIndex] = buildNodeData(scene, features, aabb, leftidx, rightidx);

        } else {
            // Current node is a leaf
            m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives);
        }
    } else {
        // Configure the current node as a giant leaf
        m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives);
    }
}

int BVH::findNextSplitLeaf()
{
    if (!m_nodes.empty()) {
        uint32_t nodecount = 0;
        std::stack<uint32_t> nodestack {};
        nodestack.push(0);

        while (!nodestack.empty()) {
            const uint32_t nodeindex = nodestack.top();
            nodestack.pop();
            const BVHInterface::Node node = m_nodes[nodeindex];

            if (node.isLeaf()) {
                if (node.primitiveCount() > BVH::LeafSize) {
                    return nodecount + 1;
                }
                nodecount++;
            } else {
                nodestack.push(node.leftChild());
                nodestack.push(node.rightChild());
            }
        }
    }
    return -1;
}

bool BVH::buildStep(const Scene& scene, const Features& features)
{
    // Find the first leaf that can be split
    uint32_t index = 0;
    bool found = false;
    if (!m_nodes.empty()) {
        std::stack<uint32_t> nodestack {};
        nodestack.push(BVH::RootIndex);

        while (!nodestack.empty()) {
            index = nodestack.top();
            const BVHInterface::Node node = m_nodes[index];
            nodestack.pop();

            if (node.isLeaf()) {
                if (node.primitiveCount() > BVH::LeafSize) {
                    found = true;
                    break;
                }
            } else {
                nodestack.push(node.leftChild());
                nodestack.push(node.rightChild());
            }
        }
    }

    if (!found)
        return false;

    const BVHInterface::Node node = m_nodes[index];
    const std::span<Primitive> prims { m_primitives.begin() + node.primitiveOffset(), node.primitiveCount() };

    // Ensure debug building
    bool old = debugBuildBVH;
    debugBuildBVH = true;

    // Convert leaf to nodes and generate 2 child leaves
    buildRecursive(scene, features, prims, index);

    debugBuildBVH = old;

    // Recalculate number of levels and leaves
    buildNumLevels();
    buildNumLeaves();

    return true;
}

// Compute the nr. of levels in your hierarchy after construction; useful for `debugDrawLevel()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLevels()
{
    uint32_t num_levels = 0;
    if (!m_nodes.empty()) {
        // BFS
        std::queue<std::tuple<uint32_t, int>> nodequeue {};
        nodequeue.emplace(BVH::RootIndex, 1);
        while (!nodequeue.empty()) {
            const uint32_t index = std::get<0>(nodequeue.front());
            const uint32_t curlevel = std::get<1>(nodequeue.front());

            const BVHInterface::Node node = m_nodes[index];
            nodequeue.pop();

            num_levels = std::max(num_levels, curlevel);

            if (!node.isLeaf()) {
                nodequeue.emplace(node.leftChild(), curlevel + 1);
                nodequeue.emplace(node.rightChild(), curlevel + 1);
            }
        }
    }
    m_numLevels = num_levels;
}

// Compute the nr. of leaves in your hierarchy after construction; useful for `debugDrawLeaf()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLeaves()
{
    uint32_t count = 0;
    if (!m_nodes.empty()) {
        std::stack<uint32_t> nodestack {};
        nodestack.push(0);

        while (!nodestack.empty()) {
            const BVHInterface::Node node = m_nodes[nodestack.top()];
            nodestack.pop();

            if (node.isLeaf()) {
                count++;
            } else {
                nodestack.push(node.leftChild());
                nodestack.push(node.rightChild());
            }
        }
    }
    m_numLeaves = count;
}

// Draw the bounding boxes of the nodes at the selected level. Use this function to visualize nodes
// for debugging. You may wish to implement `buildNumLevels()` first. We suggest drawing the AABB
// of all nodes on the selected level.
// You are free to modify this function's signature.
void BVH::debugDrawLevel(int level)
{
    // Example showing how to draw an AABB as a (white) wireframe box.
    // Hint: use draw functions (see `draw.h`) to draw the contained boxes with different
    // colors, transparencies, etc.
    if (!m_nodes.empty()) {
        // BFS
        std::queue<std::tuple<uint32_t, int>> nodequeue {};
        nodequeue.emplace(BVH::RootIndex, level);
        while (!nodequeue.empty()) {
            const uint32_t index = std::get<0>(nodequeue.front());
            const uint32_t curlevel = std::get<1>(nodequeue.front());

            const BVHInterface::Node node = m_nodes[index];
            nodequeue.pop();
            if (curlevel == 0) {
                drawAABB(node.aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.5f);
            } else if (!node.isLeaf()) {
                nodequeue.emplace(node.leftChild(), curlevel - 1);
                nodequeue.emplace(node.rightChild(), curlevel - 1);
            }
        }
    }
}

void BVH::calculateDebugSplitPlanes(std::span<BVHInterface::Primitive> nodeprims, const AxisAlignedBox& aabb)
{
    // Calculate split planes for all 3 axis
    debugSplitPlanes.clear();
    for (int i = 0; i < 3; i++) {
        std::ranges::copy(calculateSplitPlanes(aabb, i, nodeprims), std::back_inserter(debugSplitPlanes));
    }

    // Calculate costs for all split planes
    std::vector<float> costs(debugSplitPlanes.size());
    std::ranges::transform(debugSplitPlanes, costs.begin(), [](const SplitPlane& a) { return a.cost; });
    // Remove costs that have a value of inf,
    std::ranges::remove_if(costs, [](const float cost) { return !std::isnormal(cost); });

    float maxcost = *std::ranges::max_element(costs);
    float mincost = *std::ranges::min_element(costs);

    // Set color and snap split point
    for (SplitPlane& splitplane : debugSplitPlanes) {
        if (snap_splitplanes) {
            splitplane.splitpoint = computePrimitiveCentroid(nodeprims[splitplane.index]);
        }

        glm::vec4 color;
        if (std::isnormal(splitplane.cost)) {
            // Map lowest cost to 0 and highest cost to 1
            const float normalizedcost = std::clamp((splitplane.cost - mincost) / (maxcost - mincost), 0.0f, 1.0f);
            // Apply a curve to the cost to enlarge the contract in the lower cost ranges
            const float displaycost = std::pow(normalizedcost, 0.2f);
            // Map low cost to green and high cost to red
            color = { displaycost, 1.0f - displaycost, 0.05f, 0.3f };
        } else {
            // Gray if cost is inf or nan
            color = glm::vec4 { 0.5f };
        }
        splitplane.color = color;
    }
}

// Draw data of the leaf at the selected index. Use this function to visualize leaf nodes
// for debugging. You may wish to implement `buildNumLeaves()` first. We suggest drawing the AABB
// of the selected leaf, and then its underlying primitives with different colors.
// - leafIndex; index of the selected leaf.
//              (Hint: not the index of the i-th node, but of the i-th leaf!)
// You are free to modify this function's signature.
void BVH::debugDrawLeaf(int leafIndex)
{
    if (leafIndex < 0)
        return;

    const std::span<BVHInterface::Primitive> allprims = primitives();

    uint32_t nodecount = 0;
    if (!m_nodes.empty()) {
        std::stack<uint32_t> nodestack {};
        nodestack.push(0);

        while (!nodestack.empty()) {
            const uint32_t nodeindex = nodestack.top();
            nodestack.pop();
            const BVHInterface::Node node = m_nodes[nodeindex];

            if (node.isLeaf()) {
                if (nodecount == leafIndex - 1) {
                    // Draw green wireframe boundingbox
                    drawAABB(node.aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.6f);

                    const uint32_t count = node.primitiveCount();
                    const uint32_t offset = node.primitiveOffset();

                    std::span<Primitive> nodeprims(allprims.begin() + offset, count);

                    // Set values for GUI
                    num_nodes_in_primitive = count;
                    num_splitplanes = num_bins(count) - 1;

                    // Draw primitives in leaf with a half-transparent blue fill
                    glPushAttrib(GL_ALL_ATTRIB_BITS);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    glColor4f(0.05f, 0.05f, 1.0f, 0.5f);
                    for (Primitive& primitive : nodeprims) {
                        const auto& [v0, v1, v2] = std::tie(primitive.v0, primitive.v1, primitive.v2);
                        drawTriangle(v0, v1, v2);
                    }
                    glPopAttrib();

                    if (debugBuildBVH) {
                        // Calculate new split planes only if index changed
                        if (previousDebugNodeIndex != nodeindex) {
                            calculateDebugSplitPlanes(nodeprims, node.aabb);
                            previousDebugNodeIndex = nodeindex;
                        }

                        // Draw split planes
                        if (!debugSplitPlanes.empty()) {
                            for (const SplitPlane& splitplane : debugSplitPlanes) {
                                drawSplitPlane(node.aabb, splitplane.splitpoint, splitplane.axis, splitplane.color);
                            }
                        }
                    }
                }
                nodecount++;
            } else {
                nodestack.push(node.leftChild());
                nodestack.push(node.rightChild());
            }
        }
    }
}