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
    int j = image.height * texCoord.y;
    int index = (j * image.width) + i;
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
    int i1 = floor(image.width * texCoord.x);
    int i2 = ceil(image.width * texCoord.x); 
    int j1 = floor(image.height * texCoord.y);
    int j2 = ceil(image.height * texCoord.y);
    //Get average over 4
    int index1 =  (j1 * image.width) + i1;
    int index2 = (j2 * image.width) + i1;
    int index3 = (j1 * image.width) + i2;
    int index4 = (j2 * image.width) + i2;

    glm::vec3 x = (image.pixels[index1] + image.pixels[index2]) / 2.0f;
    glm::vec3 y = (image.pixels[index3] + image.pixels[index4]) / 2.0f;
        
    return (x + y) / 2.0f;
}