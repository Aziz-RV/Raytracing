#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.h"
#include "vec3.h"
#include "material.h"


struct hit_record {
    vec3 p;
    vec3 normal;
    double t;
    material* mat_ptr;
    
};

class hittable {
public:
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};

#endif
#pragma once
