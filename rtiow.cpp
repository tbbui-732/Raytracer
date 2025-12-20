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
    vec3(const f64 _x, const f64 _y, const f64 _z) : x(_x), y(_y), z(_z) {}
    vec3(const vec3* other) : x(other->x), y(other->y), z(other->z) {}

    // additional accessors for colors
    f64 r() { return x; }
    f64 g() { return y; }
    f64 b() { return z; }

    void print() { std::cout << x << ' ' << y << ' ' << z << '\n'; }

    // arithmetic
    vec3 operator+(const vec3 other) const {
        return vec3(
            this->x + other.x,
            this->y + other.y,
            this->z + other.z
        );
    }

    vec3 operator-(const vec3 other) const {
        return vec3(
            this->x - other.x,
            this->y - other.y,
            this->z - other.z
        );
    }

    vec3 operator*(const f64 scalar) const {
        return vec3(
            this->x * scalar,
            this->y * scalar,
            this->z * scalar
        );
    }

    vec3 operator/(const f64 scalar) const {
        return vec3(
            this->x / scalar,
            this->y / scalar,
            this->z / scalar
        );
    }

    void operator+=(const vec3 other) {
        x += other.x;
        y += other.y;
        z += other.z;
    }

    void operator-=(const vec3 other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
    }

    void operator*=(const f64 scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
    }

    void operator/=(const f64 scalar) {
        x /= scalar;
        y /= scalar;
        z /= scalar;
    }

    f64 length_squared() const {
        f64 x_squared = pow(this->x, 2);
        f64 y_squared = pow(this->y, 2);
        f64 z_squared = pow(this->z, 2);
        return x_squared + y_squared + z_squared;
    }

    f64 length() const {
        return sqrt(this->length_squared());
    }

    f64 magnitude() const {
        return this->length();
    }

    f64 norm() const {
        return this->length();
    }

    vec3 unit() const {
        return *this / this->length();
    }
};

vec3 operator*(const f64 scalar, const vec3 vector) {
    return vector * scalar;
}

vec3 operator/(const f64 scalar, const vec3 vector) {
    return vector / scalar;
}

f64 dot(const vec3 self, const vec3 other) {
    return (self.x * other.x) +
           (self.y * other.y) +
           (self.z * other.z);
}

vec3 cross(const vec3 self, const vec3 other) {
    return vec3(
        (self.y * other.z) - (self.z * other.y),
        (self.z * other.x) - (self.x * other.z),
        (self.x * other.y) - (self.y * other.x)
    );
}

/* END OF VEC3 */


/****************************/
/* START OF COLOR UTILITIES */
/****************************/
void write_color(const vec3 pixel_color) {
    f64 r = pixel_color.x;
    f64 g = pixel_color.y;
    f64 b = pixel_color.z;

    // translate [0,1] components to a byte range [0, 255]
    u32 rbyte = (u32)(255.999 * r);
    u32 gbyte = (u32)(255.999 * g);
    u32 bbyte = (u32)(255.999 * b);

    std::cout << rbyte << ' ' << gbyte << ' ' << bbyte << '\n';
}
/* END OF COLOR UTILITIES */

/*********************/
/* START OF RAY/RAYS */
/*********************/
struct ray {
    vec3 origin, direction;

    ray() {
        origin = vec3();
        direction = vec3();
    }

    ray(const vec3 _origin, const vec3 _direction) {
        origin = vec3(_origin);
        direction = vec3(_direction);
    }
};

/*
 * P(t) = A + t*b;
 * P(t) = Point on a line/ray
 * `A` = Origin -> typically the camera position
 * `b` = Direction
 * `t` = Some real value that dictates how far you move along `b`
 */
vec3 at(const ray r, const double t) {
    return r.origin + (r.direction * t);
}

vec3 ray_color(const ray r) {
    vec3 unit_direction = r.direction.unit();
    f64 a = 0.5 * (unit_direction.y + 1.0);
    vec3 color1 = {1.0, 1.0, 1.0};
    vec3 color2 = {0.5, 0.7, 1.0};
    return (color1 * (1.0-a)) + (color2 * a);
}
/* END OF RAY/RAYS */

int main(void) {
    // Image

    auto aspect_ratio = 16.0 / 9.0;
    int image_width = 400;

    // Calculate the image height, and ensure that it's at least 1.
    int image_height = int(image_width / aspect_ratio);
    image_height = (image_height < 1) ? 1 : image_height;

    // Camera

    auto focal_length = 1.0;
    auto viewport_height = 2.0;
    auto viewport_width = viewport_height * (double(image_width)/image_height);
    auto camera_center = vec3(0, 0, 0);

    // Calculate the vectors across the horizontal and down the vertical viewport edges.
    auto viewport_u = vec3(viewport_width, 0, 0);
    auto viewport_v = vec3(0, -viewport_height, 0);

    // Calculate the horizontal and vertical delta vectors from pixel to pixel.
    auto pixel_delta_u = viewport_u / image_width;
    auto pixel_delta_v = viewport_v / image_height;

    // Calculate the location of the upper left pixel.
    auto viewport_upper_left = camera_center
        - vec3(0, 0, focal_length) - viewport_u/2 - viewport_v/2;
    auto pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

    // Render

    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";
    for (int j = 0; j < image_height; j++) {
        for (int i = 0; i < image_width; i++) {
            auto pixel_center = pixel00_loc + (i * pixel_delta_u) + (j * pixel_delta_v);
            auto ray_direction = pixel_center - camera_center;
            ray r(camera_center, ray_direction);

            vec3 pixel_color = ray_color(r);
            write_color(pixel_color);
        }
    }
}
