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
    float i = image.width * texCoord.x;
    float j = image.height * (1 - texCoord.y);

    // Check for corners
    if (i == 0 && j == 0 || i == 0 && j == image.height || i == image.width && j == 0 || i == image.width && j == image.height) {
        if (i >= image.width) {
            i--;
        }
        if (j >= image.height) {
            j--;
        }
        float index = floor(i) + (floor(j) * image.width);
        return image.pixels[index];
    }

    //Check for edges
    if (i >= image.width - 0.5f || i <= 0.5f || j >= image.height - 0.5f || j <= 0.5f) {
        if (i >= image.width) {
            i--;
        }
        if (j >= image.height) {
            j--;
        }
        float index = floor(i) + (floor(j) * image.width);
        return image.pixels[index];
    }

    float i1 = round(i) - 0.5f;
    float i2 = round(i) + 0.5f;

    float j1 = round(j) - 0.5f;
    float j2 = round(j) + 0.5f;

    float weighti1 = abs(i - i1);
    float weighti2 = abs(i2 - i);

    //Interpolate first in x direction
    glm::vec3 pix1 = image.pixels[floor(i1) + (floor(j1) * image.width)];
    glm::vec3 pix2 = image.pixels[floor(i2) + (floor(j2) * image.width)];
    glm::vec3 interpolatedX1 = weighti2 * pix1 + weighti1 * pix2;

    glm::vec3 pix3 = image.pixels[floor(j2) * image.width + floor(i1)];
    glm::vec3 pix4 = image.pixels[floor(j2) * image.width + floor(i2)];
    glm::vec3 interpolatedX2 = weighti2 * pix3 + weighti1 * pix4;

    // Now interpolate both rows
    float distJ1 = abs(j - j1);
    float distJ2 = abs(j2 - j);

    glm::vec3 interpolated = distJ2 * interpolatedX1 + distJ1 * interpolatedX2;

    return interpolated;
}