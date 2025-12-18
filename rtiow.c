#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PROJECT_NAME "RTIOW"
#define EPSILON 1e-9

/*****************/
/* START OF VEC3 */
/*****************/
typedef struct {
    double x, y, z;
} vec3;

void print(const vec3* Vector) {
    fprintf(stdout, "%f %f %f\n", Vector->x, Vector->y, Vector->z);
}

vec3 add(const vec3* self, const vec3* other) {
    vec3 new = {
        self->x + other->x,
        self->y + other->y,
        self->z + other->z
    };
    return new;
}

vec3 sub(const vec3* self, const vec3* other) {
    vec3 new = {
        self->x - other->x,
        self->y - other->y,
        self->z - other->z
    };
    return new;
}

vec3 multiply(const vec3* Vector, const double scalar) {
    vec3 new = {
        Vector->x * scalar,
        Vector->y * scalar,
        Vector->z * scalar
    };
    return new;
}

vec3 divide(const vec3* Vector, const double scalar) {
    vec3 new = {
        Vector->x / scalar,
        Vector->y / scalar,
        Vector->z / scalar
    };
    return new;
}

double length_squared(const vec3* Vector) {
    return (Vector->x * Vector->x) +
           (Vector->y * Vector->y) +
           (Vector->z * Vector->z);
}

double length(const vec3* Vector) {
    return sqrt(length_squared(Vector));
}

double dot(const vec3* self, const vec3* other) {
    return self->x * other->x +
           self->y * other->y +
           self->z * other->z;
}

vec3 cross(const vec3* self, const vec3* other) {
    vec3 new = {
        self->y * other->z - self->z * other->y,
        self->z * other->x - self->x * other->z,
        self->x * other->y - self->y * other->x
    };
    return new;
}

vec3 unit(const vec3* Vector) {
    return divide(Vector, length(Vector));
}
/* END OF VEC3 */


/****************************/
/* START OF COLOR UTILITIES */
/****************************/
typedef vec3 color;

void write_color(const color* pixel_color) {
    double r = pixel_color->x;
    double g = pixel_color->y;
    double b = pixel_color->z;

    // translate [0,1] components to a byte range [0, 255]
    int rbyte = (int)(255.999 * r);
    int gbyte = (int)(255.999 * g);
    int bbyte = (int)(255.999 * b);

    fprintf(stdout, "%d %d %d\n", rbyte, gbyte, bbyte);
}

/* END OF COLOR UTILITIES */


int main(void) {
    const int width = 256;
    const int height = 256;

    FILE* file = fopen("image.ppm", "w");

    if (file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", "image.ppm");
        return EXIT_FAILURE;
    }

    fprintf(stdout, "P3\n%d %d\n255\n", width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            color pixel_color = {
                (double)(i) / (width-1),
                (double)(j) / (height-1),
                (double)(i) / (width-1)
            };
            write_color(&pixel_color);
        }
    }

    fclose(file);

    return EXIT_SUCCESS;
}
