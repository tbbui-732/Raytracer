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
void write_color(const vec3* pixel_color) {
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

/*********************/
/* START OF RAY/RAYS */
/*********************/
typedef struct {
    vec3 origin;
    vec3 direction;
} ray;

vec3 at(const ray* r, const double t) {
    return add(r->origin, multiply(&(r->direction), t));
}

vec3 ray_color(const ray* r) {
    vec3 new = {0, 0, 0};
    return new;
}
/* END OF RAY/RAYS */

int main(void) {
    const double aspect_ratio = 16.0 / 9.0;

    // calculate screen width and height, ensuring the height is atleast 1
    int image_width = 400;
    int image_height = (int)(image_width / aspect_ratio);
    if (image_height < 1) image_height = 1;

    // camera
    double focal_length = 1.0;
    vec3 camera_center = {0};

    // viewport
    double viewport_height = 2.0;
    double viewport_width = viewport_height *
                            (double)image_width / (double)image_height;

    // vectors across horizontal and vertical viewport edge
    vec3 viewport_u = {viewport_width, 0, 0};
    vec3 viewport_v = {0, -viewport_height, 0};

    // horizontal and vertical delta vectors from pixel to pixel
    vec3 pixel_delta_u = divide(viewport_u, width);

    // Calculate the location of the upper left pixel
    vec3 viewport_upper_left = camera_center -
                                (vec3){0, 0, focal_length} -
                                divide(viewport_u, 2) -
                                divide(viewport_u, 2);

    vec3 pixel00_location = viewport_upper_left + 0.5 *
                             add(pixel_delta_u, pixel_delta_v);

    FILE* file = fopen("image.ppm", "w");

    if (file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", "image.ppm");
        return EXIT_FAILURE;
    }

    fprintf(stdout, "P3\n%d %d\n255\n", width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            vec3 pixel_center = pixel00_location +
                                 (i * pixel_delta_u) + (j * pixel_delta_v);
            vec3 ray_direction = pixel_center - camera_center;
            ray r = {camera_center, ray_direction};
            vec3 pixel_color = ray_color(&r);
            write_color(&pixel_color);
        }
    }

    fclose(file);

    return EXIT_SUCCESS;
}
