#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cmath>

#define PROJECT_NAME "RTIOW"
#define EPSILON 1e-9

using i32 = std::int32_t;
using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using f32 = float;
using f64 = double;

// FIXME(bao): image.ppm is resulting in a white image

/*****************/
/* START OF VEC3 */
/*****************/
struct vec3 {
    f64 x, y, z;

    vec3() { x = 0.; y = 0.; z = 0.; }
    vec3(double _x, double _y, double _z) : x(_y), y(_y), z(_z) {}

    // additional accessors for colors
    double r() { return x; }
    double g() { return y; }
    double b() { return z; }

    void print() { std::cout << x << ' ' << y << ' ' << z << '\n' }

    // arithmetic
    vec3 operator+(const vec3* other) {
        return = vec3(
            this->x + other->x,
            this->y + other->y,
            this->z + other->z
        );
    }

    vec3 operator-(const vec3* other) {
        return vec3(
            this->x - other->x,
            this->y - other->y,
            this->z - other->z
        );
    }

    vec3 operator*(const double scalar) {
        return vec3(
            this->x * scalar,
            this->y * scalar,
            this->z * scalar
        );
    }

    vec3 operator/(const double scalar) {
        return vec3(
            this->x / scalar,
            this->y / scalar,
            this->z / scalar
        );
    }

    void operator+=(const vec3* other) {
        x += other->x;
        y += other->y;
        z += other->z;
    }

    void operator-=(const vec3* other) {
        x -= other->x;
        y -= other->y;
        z -= other->z;
    }

    void operator*=(const double scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
    }

    void operator/=(const double scalar) {
        x /= scalar;
        y /= scalar;
        z /= scalar;
    }

    double length_squared() {
        double x_squared = pow(this->x, 2);
        double y_squared = pow(this->y, 2);
        double z_squared = pow(this->z, 2);
        return x_squared + y_squared + z_squared;
    }

    double length() {
        return sqrt(this->length_squared());
    }

    double magnitude() {
        return this->length();
    }

    double norm() {
        return this->length();
    }

    vec3 unit() {
        return this / this->length();
    }
};

double dot(const vec3* self, const vec3* other) {
    return self->x * other->x +
           self->y * other->y +
           self->z * other->z;
}

vec3 cross(const vec3* self, const vec3* other) {
    return vec3(
        self->y * other->z - self->z * other->y,
        self->z * other->x - self->x * other->z,
        self->x * other->y - self->y * other->x
    );
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
struct ray {
    vec3 origin;
    vec3 direction;
};

vec3 at(const ray* r, const double t) {
    /*
        P(t) = A + t*b;
        s.t. `A` is the origin, `b` is the direction, and `t` is real-value
    */
    vec3 t_times_b = multiply(&r->direction, t);
    return add(&r->origin, &t_times_b);
}

vec3 ray_color(const ray* r) {
    vec3 unit_direction = unit(&r->direction);
    double a = 0.5 * (unit_direction.y + 1.0);

    vec3 color1 = {1.0, 1.0, 1.0};
    vec3 color2 = {0.5, 0.7, 1.0};

    vec3 op1 = multiply(&color1, (1.0-a));
    vec3 op2 = multiply(&color2, a);
    vec3 result = add(&op1, &op2);
    return result;
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
    vec3 pixel_delta_u = divide(&viewport_u, (double)image_width);
    vec3 pixel_delta_v = divide(&viewport_v, (double)image_height);

    // Calculate the location of the upper left pixel
    vec3 viewport_u_half = divide(&viewport_u, 2);
    vec3 viewport_v_half = divide(&viewport_v, 2);
    vec3 viewport_dif = sub(&viewport_u_half, &viewport_v_half);

    vec3 viewport_upper_left = sub(&camera_center, &(vec3){0, 0, focal_length});
    viewport_upper_left = sub(&viewport_upper_left, &viewport_dif);

    vec3 pixel00_location = add(&pixel_delta_u, &pixel_delta_v);
    pixel00_location = multiply(&pixel00_location, 0.5);
    pixel00_location = add(&viewport_upper_left, &pixel00_location);

    if (file == nullptr) {
        std::cerr << "Error opening files: image.ppm\n";
        return EXIT_FAILURE;
    }

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n"
    for (int j = 0; j < image_height; j++) {
        for (int i = 0; i < image_width; i++) {
            vec3 iu = multiply(&pixel_delta_u, i);
            vec3 jv = multiply(&pixel_delta_v, j);
            vec3 pixel_center = add(&iu, &jv);
            pixel_center = add(&pixel_center, &pixel00_location);

            vec3 ray_direction = sub(&pixel_center, &camera_center);
            ray r = {camera_center, ray_direction};
            vec3 pixel_color = ray_color(&r);
            write_color(&pixel_color);
        }
    }

    return EXIT_SUCCESS;
}
