#include "texture.h"
#include "render.h"
#include <framework/image.h>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    int i = image.width * texCoord.x;
    int j = image.height * (1 - texCoord.y);

    if (i == image.width) {
        i--;
    }
    if (j == image.height) {
        j--;
    }

    float index = floor(j) * image.width + floor(i);
    return image.pixels[index];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    float i = (image.width * texCoord.x);
    float j = (image.height * texCoord.y);

    //Get 4 courners
    float x1 = (round(i) - 0.5f) / image.width;
    float y1 = (round(j) - 0.5f) / image.height;
    float x2 = (round(i) + 0.5f) / image.width;
    float y2 = (round(j) + 0.5f) / image.height;

    float weight1 = abs(texCoord.x - x1) * image.width;
    float weight2 = abs(texCoord.y - y1) * image.height;

    //Get closest 
     glm::vec3 a1 = sampleTextureNearest(image, { x1, y1 });
     glm::vec3 a2 = sampleTextureNearest(image, { x2, y1 });
     glm::vec3 b1 = sampleTextureNearest(image, { x1, y2 });
     glm::vec3 b2 = sampleTextureNearest(image, { x2, y2 });


    //interpolate in x direction
    glm::vec3 interpolatedx1 = a1 * (1 - weight1) + (weight1 * a2);
    glm::vec3 interpolatedx2 = (b2 * weight1) + (1 - weight1) * b1; 

    //Interpolate in y direction
    glm::vec3 interpolated = weight2 * interpolatedx2 + (1 - weight2) * interpolatedx1;
    return interpolated;

}
