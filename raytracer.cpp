#include <iostream>
#include <fstream>
#include <cmath>

struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator/(double t) const { return Vec3(x / t, y / t, z / t); }

    double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3 normalize() const { double len = std::sqrt(x * x + y * y + z * z); return *this / len; }
};

struct Ray {
    Vec3 origin, direction;

    Ray(const Vec3& origin, const Vec3& direction) : origin(origin), direction(direction.normalize()) {}
};

bool hit_sphere(const Vec3& center, double radius, const Ray& ray, double& t) {
    Vec3 oc = ray.origin - center;
    double a = ray.direction.dot(ray.direction);
    double b = 2.0 * oc.dot(ray.direction);
    double c = oc.dot(oc) - radius * radius;
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return false;
    t = (-b - std::sqrt(discriminant)) / (2.0 * a);
    return true;
}

bool hit_plane(const Vec3& point, const Vec3& normal, const Ray& ray, double& t) {
    double denom = normal.dot(ray.direction);
    if (std::fabs(denom) > 1e-6) { // Avoid division by zero
        Vec3 p0l0 = point - ray.origin;
        t = p0l0.dot(normal) / denom;
        return t >= 0;
    }
    return false;
}

Vec3 ray_color(const Ray& ray) {
    Vec3 sphere_center(0, 0, -1);
    double sphere_radius = 0.5;
    double t;
    if (hit_sphere(sphere_center, sphere_radius, ray, t)) {
        return Vec3(1, 0, 0); // Red sphere
    }

    Vec3 plane_point(0, -0.5, -1);
    Vec3 plane_normal(0, 1, 0);
    if (hit_plane(plane_point, plane_normal, ray, t)) {
        return Vec3(0.5, 0.5, 0.5); // Gray plane
    }

    // Background gradient
    Vec3 unit_direction = ray.direction.normalize();
    t = 0.5 * (unit_direction.y + 1.0);
    return Vec3(1.0, 1.0, 1.0) * (1.0 - t) + Vec3(0.5, 0.7, 1.0) * t;
}

void render_scene(const std::string& filename, int image_width, int image_height) {
    std::ofstream file(filename);
    file << "P3\n" << image_width << " " << image_height << "\n255\n";

    Vec3 origin(0, 0, 0);
    double viewport_height = 2.0;
    double viewport_width = 2.0 * image_width / image_height;
    double focal_length = 1.0;

    Vec3 horizontal(viewport_width, 0, 0);
    Vec3 vertical(0, viewport_height, 0);
    Vec3 lower_left_corner = origin - horizontal / 2 - vertical / 2 - Vec3(0, 0, focal_length);

    for (int j = image_height - 1; j >= 0; --j) {
        for (int i = 0; i < image_width; ++i) {
            double u = double(i) / (image_width - 1);
            double v = double(j) / (image_height - 1);
            Ray ray(origin, lower_left_corner + horizontal * u + vertical * v - origin);
            Vec3 color = ray_color(ray);
            file << static_cast<int>(255.999 * color.x) << " "
                << static_cast<int>(255.999 * color.y) << " "
                << static_cast<int>(255.999 * color.z) << "\n";
        }
    }
    file.close();
}

int main() {
    const int image_width = 400;
    const int image_height = 200;
    render_scene("output.ppm", image_width, image_height);
    return 0;
}