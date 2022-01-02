#include <Windows.h>
#include <iostream>
#include <fstream>
#include <ostream>
#include <random>
#include <stdlib.h>
#include <stdio.h>
#include "vec3.h"
#include "color.h"
#include "ray.h"
#include "hittable.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "float.h"
#include "material.h"
using namespace std;

double hit_sphere(const vec3& center, double radius, const ray& r) {

    vec3 oc = r.origin() - center;
    double a = dot(r.direction(), r.direction());
    double b = 2.0 * dot(oc, r.direction());
    double c = dot(oc, oc) - radius * radius;
    double discriminant = b * b - 4 * a * c;


    if (discriminant < 0) {
        return -1.0;
    }
    else {
        return (-b - sqrt(discriminant)) / (2.0 * a);
    }

}
bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;

    for (int i = 0; i < list_size; i++) {
        if (list[i]->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}
bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) {
    vec3 oc = r.origin() - center;
    double a = dot(r.direction(), r.direction());
    double b = dot(oc, r.direction());
    double c = dot(oc, oc) - radius * radius;
    double discriminant = b * b - a * c;
    if (discriminant > 0) {
        double temp = (-b - sqrt(b * b - a * c)) / a;
        if (temp < t_max && t_min < temp) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;



            return true;
        }


        // Find the nearest root that lies in the acceptable range.
        temp = (-b - sqrt(b * b - a * c)) / a;
        if (temp < t_max && t_min < temp) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;

            return true;
        }



    }
    return false;
}


vec3 ray_color(const ray& r) {

    double t = hit_sphere(vec3(0, 0, -1), 0.5, r);
    if (t > 0.0) {
        vec3 N = unit_vector(r.at(t) - vec3(0, 0, -1));
        return 0.5 * vec3(N.x() + 1, N.y() + 1, N.z() + 1);
    }

    vec3 unit_direction = unit_vector(r.direction());

    t = 0.5 * (unit_direction.y() + 1.0);

    return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
}
vec3 ray_color(const ray& r, hittable* world) {
    hit_record rec;
    if (world->hit(r, 0.001, RAND_MAX, rec))
    {
        vec3 target = rec.p + rec.normal + random_in_unit_sphere();
        return 0.5 * ray_color(ray(rec.p, target - rec.p), world);
    }
    else
    {
        vec3 unit_direction = unit_vector(r.direction());
        auto t = 0.5 * (unit_direction.y() + 1.0);

        return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
    }


}
vec3 ray_color(const ray& r, hittable* world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.

    if (world->hit(r, 0.001, RAND_MAX, rec)) {

        ray scattered;
        vec3 attenuation;
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            return attenuation * ray_color(scattered, world, depth + 1);
        }

        else {
            return vec3(0, 0, 0);
        }

    }
    else {

        vec3 unit_direction = unit_vector(r.direction());
        auto t = 0.5 * (unit_direction.y() + 1.0);
        return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
    }
}
hittable *random_scene() {
    int n = 500;
    sphere **list = new sphere*[n + 1];
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));
    int i = 1;

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat =(double)rand() / RAND_MAX;
            point3 center(a + 0.9 * (rand() / RAND_MAX), 0.2, b + 0.9 * (rand() / RAND_MAX));

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
              

                if (choose_mat < 0.8) {
                    // diffuse
                    list[i++] = new sphere(center, 0.2, new lambertian(vec3((rand() / RAND_MAX) * (rand() / RAND_MAX), (rand() / RAND_MAX) * (rand() / RAND_MAX), (rand() / RAND_MAX) * (rand() / RAND_MAX))));
                        
                }
                else if (choose_mat < 0.95) {
                    // metal
                  
                    list[i++] = new sphere(center, 0.2, new metal(vec3(0.5*(1+ (rand() / RAND_MAX)), 0.5 * (1 + (rand() / RAND_MAX)), 0.5 * (1 + (rand() / RAND_MAX))),0.5*(rand() / RAND_MAX)));
                }
                else {
                    // glass
                    list[i++] = new sphere(center, 0.2, new dielectric(1.5));
                }
            }
        }
    }

   list[i++]=new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));

    list[i++]=new sphere(vec3(-4, 1, 0), 1.0,new lambertian (vec3(0.4, 0.2, 0.1)));

    list[i++]=new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

    return new hittable_list(list,i);
}


int main() {

    // Image
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int ns = 100;

    // Camera

    auto viewport_height = 2.0;
    auto viewport_width = aspect_ratio * viewport_height;
    auto focal_length = 1.0;

    auto origin = vec3(0, 0, 0);
    auto horizontal = vec3(viewport_width, 0, 0);
    auto vertical = vec3(0, viewport_height, 0);
    auto lower_left_corner = origin - horizontal / 2 - vertical / 2 - vec3(0, 0, focal_length);

    ofstream outputfile;
    outputfile.open("imageRaytracing.ppm");

    // Render

    outputfile << "P3\n" << image_width << ' ' << image_height << "\n255\n";
  
    auto world = random_scene();

    vec3 lookfrom(13, 3, 2);
    vec3 lookat(0, 0, -1);
    double dist_to_focus=(lookfrom - lookat).length();
    double aperture = 2.0;
    
    camera cam(lookfrom,lookat,vec3(0,1,0),20,double(image_height)/double(image_width),aperture,dist_to_focus);


 
    for (int j = image_height - 1; j >= 0; j--) {
       for(int i = 0; i < image_width; i++) {
            vec3 col(0, 0, 0);
            for (int s = 0; s < ns; s++) {
                // mon projet est sous windows donc j'ai remplace drand48 par rand() / RAND_MAX
                double ru = (double)rand() / RAND_MAX;
                double rv = (double)rand() / RAND_MAX;
                double u = double(i + ru) / double(image_width);
                double v = double(j + rv) / double(image_height);
                ray r = cam.get_ray(u, v);
                vec3 p = r.at(2.0);
                col += ray_color(r, world,0);

            }
            col /= double(ns);
            col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
            int ir = static_cast<int>(255.999 * col[0]);
            int ig = static_cast<int>(255.999 * col[1]);
            int ib = static_cast<int>(255.999 * col[2]);

            outputfile << ir << ' ' << ig << ' ' << ib << '\n';
        }
    }
    outputfile.close();
}