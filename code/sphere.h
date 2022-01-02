#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"
#include <cmath>
#include "material.h"


class sphere : public hittable {
public:
    sphere() {}
    sphere(vec3 cen, double r)
        : center(cen), radius(r){};
    sphere(vec3 cen, double r,material *m)
        : center(cen), radius(r),mat_ptr(m) {};

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const ;
    vec3 center;
    double radius;
    material *mat_ptr;
    
};

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
#endif
