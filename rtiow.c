#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PROJECT_NAME "RTIOW"
#define EPSILON 1e-9

/* START OF VEC3 */
typedef struct {
    double x, y, z;
} vec3;

void print(const vec3* Vector) {
    printf("%f %f %f\n", Vector->x, Vector->y, Vector->z);
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

int main(void) {
    const int width = 256;
    const int height = 256;

    FILE* file = fopen("image.ppm", "w");

    if (file == NULL) {
        printf("Error opening file: %s\n", "image.ppm");
        return EXIT_FAILURE;
    }

    fprintf(file, "P3\n%d %d\n255\n", width, height);

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            double r = (double)(i) / (width-1);
            double g = (double)(j) / (height - 1);
            double b = 0.0;

            int ir = (int)(255.999 * r);
            int ig = (int)(255.999 * g);
            int ib = (int)(255.999 * b);

            fprintf(file, "%d %d %d\n", ir, ig, ib);
        }
    }

    fclose(file);

    return EXIT_SUCCESS;
}
