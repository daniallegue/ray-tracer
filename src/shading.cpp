#include "render.h"
#include "texture.h"
#include <cmath>
#include <fmt/core.h>
#include <glm/geometric.hpp>
#include <glm/gtx/string_cast.hpp>
#include <shading.h>

// This function is provided as-is. You do not have to implement it (unless
// you need to for some extra feature).
// Given render state and an intersection, based on render settings, sample
// the underlying material data in the expected manner.
glm::vec3 sampleMaterialKd(RenderState& state, const HitInfo& hitInfo)
{
    if (state.features.enableTextureMapping && hitInfo.material.kdTexture) {
        if (state.features.enableBilinearTextureFiltering) {
            return sampleTextureBilinear(*hitInfo.material.kdTexture, hitInfo.texCoord);
        } else {
            return sampleTextureNearest(*hitInfo.material.kdTexture, hitInfo.texCoord);
        }
    } else {
        return hitInfo.material.kd;
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the scene-selected shading model, returning the reflected light towards the target.
glm::vec3 computeShading(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Hardcoded linear gradient. Feel free to modify this
    static LinearGradient gradient = {
        .components = {
            { 0.1f, glm::vec3(215.f / 256.f, 210.f / 256.f, 203.f / 256.f) },
            { 0.22f, glm::vec3(250.f / 256.f, 250.f / 256.f, 240.f / 256.f) },
            { 0.5f, glm::vec3(145.f / 256.f, 170.f / 256.f, 175.f / 256.f) },
            { 0.78f, glm::vec3(255.f / 256.f, 250.f / 256.f, 205.f / 256.f) },
            { 0.9f, glm::vec3(170.f / 256.f, 170.f / 256.f, 170.f / 256.f) },
        }
    };

    if (state.features.enableShading) {
        switch (state.features.shadingModel) {
            case ShadingModel::Lambertian:
                return computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::Phong:
                return computePhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::BlinnPhong:
                return computeBlinnPhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::LinearGradient:
                return computeLinearGradientModel(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
        };
    }

    return lightColor * sampleMaterialKd(state, hitInfo);
}

// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a Lambertian diffuse shading, returning the reflected light towards the target.
glm::vec3 computeLambertianModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    float angle = glm::dot(glm::normalize(hitInfo.normal), glm::normalize(lightDirection));
    if (angle <= 0.0) {
        return glm::vec3(0, 0, 0);
    };
    return lightColor * sampleMaterialKd(state, hitInfo) * angle;
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    if (glm::dot(lightDirection, hitInfo.normal) < 0) {
        return glm::vec3 { 0 };
    }

    float angle = glm::dot(glm::normalize(hitInfo.normal), glm::normalize(lightDirection));
    

    float k = 2.0f * (glm::dot((lightDirection) / glm::length(lightDirection), hitInfo.normal / glm::length(hitInfo.normal)));
    glm::vec3 r = (lightDirection) - k * hitInfo.normal;
    glm::vec3 v = (cameraDirection) / glm::length(cameraDirection);

    float product = glm::dot(r / glm::length(r), v / glm::length(v));

    return lightColor * ( hitInfo.material.ks * pow(product, hitInfo.material.shininess) / (glm::length(v) * glm::length(r)));
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Blinn-Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBlinnPhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    glm::vec3 v1 = (cameraDirection) / glm::length(cameraDirection);
    glm::vec3 v2 = (lightDirection) / glm::length(lightDirection);
    //Halfway vector
    glm::vec3 v = (v1 + v2) / glm::length(v1 + v2);
    if (glm::dot(hitInfo.normal / glm::length(hitInfo.normal), lightDirection) < 0) {
        return glm::vec3 { 0 };
    }

    return lightColor * (hitInfo.material.ks * pow(glm::dot(glm::normalize(hitInfo.normal), v), hitInfo.material.shininess));
    
}

// TODO: Standard feature
// Given a number ti between [-1, 1], sample from the gradient's components and return the
// linearly interpolated color, for which ti lies in the interval between the t-values of two
// components, or on a boundary. If ti falls outside the gradient's smallest/largest components,
// the nearest component must be sampled.
// - ti; a number between [-1, 1]
// This method is unit-tested, so do not change the function signature.
glm::vec3 LinearGradient::sample(float ti) const
{
    //Get min-max
    std::vector<float> ts;
    for (int i = 0; i < components.size(); i++) {
        ts.push_back(components[i].t);
    }
    float min = ts[0];
    float minIndex = 0;
    float max = ts[0];
    float maxIndex = 0;
    for (int i = 1; i < components.size(); i++) {
        if (components[i].t > max) {
                max = components[i].t;
                maxIndex = i;
        }
        if (components[i].t < min) {
                min = components[i].t;
                minIndex = i;
        }
    }

    //Check for out of boundaries
    if (ti > max) {
        return components[maxIndex].color;
    }
    if (ti < min) {
        return components[minIndex].color;
    }


    for (int i = 0; i < components.size() - 1; i++) {
        float t = components[i].t;
        float t2 = components[i + 1].t;
        if (t == ti) {
            return components[i].color;
              
        }
        if (t2 == ti) {
            return components[i + 1].color;
           
        }
        if (ti > t && ti < t2) {
            glm::vec3 res = { components[i].color.x + components[i + 1].color.x,
                components[i].color.y + components[i + 1].color.y,
                components[i].color.z + components[i + 1].color.z };
            return res * 0.5f; 
        }
    }



}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a diffuse shading model, such that the diffuse component is sampled not
// from the intersected material, but a provided linear gradient, based on the cosine of theta
// as defined in the diffuse shading part of the Phong model.
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - gradient;        the linear gradient object
// - return;          the result of shading
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeLinearGradientModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    float cos_theta = glm::dot(glm::normalize(lightDirection), glm::normalize(hitInfo.normal));
    glm::vec3 t = gradient.sample(cos_theta);
    return glm::vec3 { t.x * lightColor.x, t.y * lightColor.y, t.z * lightColor.z } * cos_theta;
 
}