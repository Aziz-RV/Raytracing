#ifndef MATERIAL_H
#define MATERIAL_H
#include "vec3.h"
#include "ray.h"
#include "sphere.h"
#include "hittable_list.h"
#include "hittable.h"

class material {
public:
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const = 0;
};
class lambertian : public material {
public:
    lambertian(const vec3& a) : albedo(a) {}

    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
        vec3 target = rec.normal + random_in_unit_sphere();
        scattered = ray(rec.p, target - rec.p);
        attenuation = albedo;
        return true;
    }


    vec3 albedo;
};
class metal : public material {
public:
    metal(const vec3& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}
    virtual bool scatter( const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const  {
        vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere());
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);
    }

    vec3 albedo;
    double fuzz;
};
class dielectric : public material {
public:
    dielectric(double index_of_refraction) : ir(index_of_refraction) {}

    virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const {

        vec3 outward_normal;
        vec3 reflected = reflect(r_in.direction(), rec.normal);
        double ni_over_nt;
        attenuation = vec3(1.0, 1.0, 1.0);
        vec3 refracted;
        fuble reflect_prob;
        double cosine;
        if (dot(r_in.direction(), rec.normal) > 0) {
            outward_normal = -rec.normal:
            ni_over_nt = ir;
            cosine = ir * dot(r_in.direction(), rec.normal) / r_in.direction().length();
        }
        else {
            outward_normal = rec.normal:
            ni_over_nt = 1.0 / ir;
            cosine = -dot(r_in.direction(), rec.normal) / r_in.direction().length();

        }
        if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted)) {
            reflect_prob = schlick(cosine, ir);
        }
        else {
            scattered = ray(rec.p, reflected);
            reflect_prob = 1.0;

        }
        if (((double)rand() / RAND_MAX) < reflect_prob) {
            scattered = ray(rec.p, reflected);
        }
        else {
            scattered = ray(rec.p, refracted);
        }
        return true;
    }
    double ir; // Index of Refraction
};

#endif