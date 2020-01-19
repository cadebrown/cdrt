/* render.c - the default rendering engine for raytracing 

At a basic level, a lot of ray intersection tests are done on the entire geometry of the scene

These will each result in pixel values. A higher level, stack-based recursion algorithm is used to
combine individual rays, and recursively create more rays for things that reflect/refract.

TODO: add global illumination & direct illumination/lighting





*/

#include "cdrt.h"

// minimal accuracy, that everything should be accurate to
#define E0 (1e-12f)
// slightly more tolerant error constant
#define E1 (1e-10f)
// amount to back off a ray off the surface to stop self intersection
#define E2 (1e-6f)
// amount to back off refraction rays
#define E3 (1e-2f)


// infinity
#define INF INFINITY

// pi constant
#define PI 3.141592653589793238f


/* struct RenderItem - represents a rendering job to complete 

The resulting color of the RenderItem should be (compute(ray) * mask)

*/

struct RenderItem {
    // depth in the tree of rays cast
    int depth;

    // the ray to cast out and attempt to render
    ray_t ray;

    // the mask of the color (i.e. what to multiply by)
    v3 mask;

};


/* image sampling functions */

// default image sampler given uv coordinates (0,0)->top left, (1,1)->bottom right
static inline v3 img_sample(img_t img, float u, float v) {
    // get size params
    int w = img->w, h = img->h;

    // get the pixel values as floats
    float pxf = u, pyf = (1.0f - v);

    // convert to (0, 1) range
    pxf = pxf - floorf(pxf);
    pyf = pyf - floorf(pyf);
    if (pxf < 0) pxf++;
    if (pyf < 0) pyf++;

    // convert to image scaled
    pxf *= w;
    pyf *= h;

    // integer indexes
    int px = (int)pxf, py = (int)pyf;
    // next pixel (with wraparound computed)
    int npx = px == w - 1 ? 0 : px, npy = py == h - 1 ? 0 : py;

    /* linear interpolation: take 4 samples, average them */

    // get a group of 2x2 pixels
    pix_t   tlp = img->data[py * w + px],  trp = img->data[py * w + npx],
            blp = img->data[npy * w + px], brp = img->data[npy * w + npx];

    // convert to v3
    v3  tl = V3(tlp.r, tlp.g, tlp.b), tr = V3(trp.r, trp.g, trp.b),
        bl = V3(blp.r, blp.g, blp.b), br = V3(brp.r, brp.g, brp.b);


    // now, combine rows into single values on the left side
    tl = v3_lerp(tl, tr, pxf - px);
    bl = v3_lerp(bl, br, pxf - px);

    // combine rows 
    tl = v3_lerp(tl, bl, pyf - py);

    // now, transform into (0, 1) range
    return v3_scale(tl, 1.0f / 255.f);
}


/* generic collision functions (i.e. no graphics code neccessarily) */

// return the hit results, searching up to max_dist (use INF for infinite distance)
// to check whether or not the raycast hit, check `result.hit`
static rch_t raycast(scene_t scene, ray_t ray, float max_dist) {

    // keep track of the closest hit, starting at very large values
    rch_t closest = {NULL, max_dist};

    // cached values

    // one of MODEL_* values
    int hit_type;
    // the hit index into the specific array
    int hit_idx = -1;

    // temporary variables
    bool is_inside = false;

    // iteration variable
    int i;

    // first, check if it hit any spheres
    for (i = 0; i < scene->cache.n_spheres; ++i) {
        // current sphere
        v4 cs = scene->cache.spheres[i];

        // do some geometry to test intersection, get the displacement
        v3 L = v3_sub(cs.xyz, ray.origin);

        // calculate the radius squared
        float r2 = cs.w * cs.w;
        // cos(A), where A is angle between the ray and the line connecting
        // the ray and the middle of the sphere
        float tca = v3_dot(L, ray.dir);
        // check the horizontal distance between them
        float d2 = v3_mag2(L) - tca * tca;

        // if its larger than the sphere, it didn't hit
        if (d2 > r2) continue;
        //if (r2 - d2 > tca * tca) return RCH_MISS;
        // calculate the actual distance
        // first, compute the radical R, where the d=avg(+/-)R
        float tch = sqrtf(r2 - d2);
        // solution (distances) = tca +- tch, so try both
        float dist_away = tca - tch;
        if ((is_inside = (dist_away <= E1))) dist_away = tca + tch;

        // either didn't hit, or hit too far away
        if (dist_away <= E0 || dist_away > closest.dist) continue;


        hit_idx = i;
        hit_type = MODEL_SPHERE;

        // set the distance
        closest.dist = dist_away;

        // capture this variable
        closest.is_inside = is_inside;

        // compute the rest in the later block

    }

    // next, check for planes
    for (i = 0; i < scene->cache.n_planes; ++i) {
        // current plane
        v4 cp = scene->cache.planes[i];
        // dot it with
        float dist_off = v4_dot(cp, V4(ray.origin.x, ray.origin.y, ray.origin.z, 1.0f));
        if (fabsf(dist_off) < E0) continue;

        // compute the dot product between the normal and the ray direction
        float dp = v3_dot(cp.xyz, ray.dir);

        // now calculate distance away from the view
        float dist_away = dist_off / dp;

        // check and make sure we still care about it
        if (dist_away < E0 || dist_away >= closest.dist) continue;

        hit_idx = i;
        hit_type = MODEL_PLANE;

        // update the current hit
        closest.dist = dist_away;

        // they are facing the same way
        closest.is_inside = dp > 0.0f;

        // compute the rest in the later block
    }


    if (hit_idx >= 0) {

        // we actually hit something, so fill relevant details
        // this saves computation time, since we don't compute UV/normal/etc
        // for every possible hit, but rather only the closest
        if (hit_type == MODEL_SPHERE) {

            // get the actual hit object
            obj_t hit_obj = scene->cache.ext_spheres[hit_idx].src;
            // and its data
            v4 cs = scene->cache.spheres[hit_idx];

            // set it on the rch
            closest.hit = hit_obj;

            // calculate the point at which it hit
            closest.point = v3_add(ray.origin, v3_scale(ray.dir, closest.dist));

            // calculate the normal of the surface (normalized)
            closest.normal = v3_unit(v3_sub(closest.point, cs.xyz));

            // if inside, the normal is pointing inwards as well
            closest.normal = is_inside ? v3_neg(closest.normal) : closest.normal;

            // basically a inverse mercator projection with rotation, to get the UV
            v3 diff = m3x3_mul_v3(scene->cache.ext_spheres[hit_idx].UVR, closest.normal);
            //printf("::"  M3X3_FMT "\n", M3X3__(scene->cache.ext_spheres[hit_idx].UVR));


            //v3 diff = closest.normal;
            closest.uv = V2(0.5f + atan2f(diff.z, diff.x) / (2.0f * PI), 0.5f + asinf(diff.y) / PI);
            
        } else if (hit_type == MODEL_PLANE) {

            // get the actual hit object
            obj_t hit_obj = scene->cache.ext_planes[hit_idx].src;
            // and its data
            v4 cp = scene->cache.planes[hit_idx];

            // set it on the rch
            closest.hit = hit_obj;

            closest.normal = closest.is_inside ? v3_neg(cp.xyz) : cp.xyz;
            closest.point = v3_add(ray.origin, v3_scale(ray.dir, closest.dist));

            //closest.uv = V2_(closest.point.x, closest.point.z);
            v3 localH = m3x3_mul_v3(scene->cache.ext_planes[hit_idx].UVT, closest.point);
            closest.uv = V2(localH.x, localH.z);
            //closest.uv = V2(closest.point.x, closest.point.z);
            //closest.midx = closest_plane.mat;

        } else {
            // this shouldn't happen, but don't do anythingin this case
        }
    }

    return closest;
}


/* graphics functionality */


// calculate DI (direct illumination) for a given hit
// DI is essentially the light from light sources directly, not counting
// any sub bounces
static inline v3 render_DI(scene_t scene, rch_t hit) {
    // the sum of all the colors
    v3 color = V3_0;

    // calculate directional lighting
    int i;
    for (i = 0; i < scene->cache.n_dirlights; ++i) {
        // get the current directional light
        struct dirlight_s dirl = scene->cache.dirlights[i];

        // construct a ray pointing to the light
        ray_t toL = (ray_t){ v3_add(hit.point, v3_scale(hit.normal, E3)), v3_neg(dirl.dir) };

        // check light visibility
        if (!raycast(scene, toL, INF).hit) {
            //create a dot product to see how much light
            float dp = -v3_dot(dirl.dir, hit.normal);
            if (dp > 0.0f) {
                // we are facing the light, so add a scaled version
                color = v3_add(color, v3_scale(dirl.color, dp));
            }
        }
    }

    // calculate point lighting
    for (i = 0; i < scene->cache.n_pointlights; ++i) {
        // get current point light
        struct pointlight_cache_s ptl = scene->cache.pointlights[i];

        // construct a ray pointing to the light
        ray_t toL = (ray_t){ v3_add(hit.point, v3_scale(hit.normal, E3)), v3_unit(v3_sub(ptl.pos, hit.point)) };
        
        // compute distance to light
        float dist = v3_dist(hit.point, ptl.pos);

        // check light visibility
        if (!raycast(scene, toL, dist - E2).hit) {
            // nothing was in the way, so add lighting calculations
            // use modified inverse square law, with 1/(a+d^2), where a
            // is a dampening factor
            float amt = v3_dot(toL.dir, hit.normal) / (1.0f + dist*dist);

            // add to the result
            if (amt > 0) color = v3_add(color, v3_scale(ptl.ptl.color, amt));
        }

    }

    return color;
}


// calculate GI (global illumination) for a given hit
static inline v3 render_GI(scene_t scene, rch_t hit, int samples) {
    // sum of the color
    v3 color = V3_0;

    // total sum of the weights of the measurements
    float coef = 0.0f;

    // the ray to shoot
    ray_t ray;

    // current result
    rch_t rch;

    // get square root of the number of samples
    int n = (int)(sqrtf(samples));

    int i;
    // take samples from the hemisphere
    for (i = 0; i < n; ++i) {
        int j;
        for (j = 0; j < n; ++j) {

            // sweep through
            ray.dir = v3_sphere_sample_ZA(2.0f * i / n - 1.0f, (2.0f * PI) * j / n);

            float dp = v3_dot(ray.dir, hit.normal);

            // if dp < 0 , then it is facing against the normal, so we flip it,
            // since we are sampling on the positive hemisphere
            if (dp < 0) {
                dp = -dp;
                ray.dir = v3_neg(ray.dir);
            }
            coef += dp;

            // take the origin off a bit to avoid self intersection
            ray.origin = v3_add(hit.point, v3_scale(ray.dir, E2));

            // send GI ray output
            if (!(rch = raycast(scene, ray, INF)).hit) {
                // GI bounce missed, return sky color
                color = v3_add(color, v3_scale(V3(.1, .1, .1), dp));
                continue;
            }

            // else, we hit something, so get its material
            mat_t mat = rch.hit->model.mat;

            if (mat->tex_diff) {
                // add diffuse color

                // sample the image
                v3 sample = img_sample(mat->tex_diff, rch.uv.x, rch.uv.y);

                // calculate direct illumination
                v3 di_calc = render_DI(scene, rch);

                // add it to the result
                color = v3_add(color, v3_scale(v3_mul(di_calc, sample), dp));

            }
        }
    }

    // scale the result
    return v3_scale(color, 1.0f / (PI * coef));
}



// render a single ray and return a color
static v3 render_ray(scene_t scene, ray_t primary_ray) {

    // current hit
    rch_t rch;

    // current material
    mat_t mat;

    // the accumulated color
    v3 vsum = V3_0;


    /* Essentially: keep a stack of rendering operations

    Even though the constant '256' is there

    */

    struct RenderItem rstk[256];
    
    // current index to the render stack, start at 0,
    // which is the primary ray
    int rstk_i = 0;
    rstk[0].depth = 0;

    // the main ray
    rstk[0].ray = primary_ray;
    // full color
    rstk[0].mask = V3_1;

    // number of rays processed
    int ct = 0;


    while (rstk_i >= 0 && ct < 100) {
        // pop off an item from the render stack
        struct RenderItem cur = rstk[rstk_i--];

        // attempt a raycast
        if (!(rch = raycast(scene, cur.ray, INF)).hit) {
            // add background, since we failed to hit it
            vsum = v3_add(vsum, v3_mul(cur.mask, V3(0.1f, 0.1f, 0.1f)));
            continue;
        }

        // assume we have a valid material
        if (!rch.hit->model.mat) continue;

        // else, we did hit something to render

        // get the material
        mat = rch.hit->model.mat;

        if (mat->tex_diff) {
            // sample the image at given coordinates
            v3 sample = img_sample(mat->tex_diff, rch.uv.x, rch.uv.y);

            if (cur.depth < 1) {
                // compute global illumination
                v3 GI = render_GI(scene, rch, 100);

               // printf("GI: " V3_FMT "\n", V3__(GI));

                vsum = v3_add(vsum, v3_mul(v3_mul(cur.mask, sample), GI));

            }



            // compute direct lighting
            v3 DI = render_DI(scene, rch);

            // add it to the value of the pixel
            vsum = v3_add(vsum, v3_mul(v3_mul(cur.mask, sample), DI));



        } else if (cur.depth < 8) {

            if (mat->tex_spec && mat->tex_refr) {
                // reflective and refractive material

                // refractive material
                v3 sample = img_sample(mat->tex_refr, rch.uv.x, rch.uv.y);

                // create a reflective ray
                ray_t ray_refl = (ray_t){ 
                    v3_add(rch.point, v3_scale(rch.normal, E2)), 
                    v3_refl(cur.ray.dir, rch.normal) 
                };

                // create a refraction ray
                ray_t ray_refr = (ray_t){ 
                    v3_add(rch.point, v3_scale(rch.normal, rch.is_inside ? -E2 : E2)), 
                    v3_neg(v3_refract(cur.ray.dir, rch.normal, mat->IOR))
                };

                // amount of specular color
                // default to 1 in case of total internal refraction
                float spec_mix = 1.0f;

                if (!v3_eqe(ray_refr.dir, V3_0, E2)) {
                //if (rch.is_inside) printf("is_inside: %d\n", (int)rch.is_inside);

                    // compute the specular mix
                    spec_mix = v3_fresnel(cur.ray.dir, rch.normal, mat->IOR);

                    rstk[++rstk_i] = (struct RenderItem){ cur.depth+1, ray_refr, v3_scale(sample, 1.0f - spec_mix) };

                }
                rstk[++rstk_i] = (struct RenderItem){ cur.depth+1, ray_refl, v3_scale(sample, spec_mix) };

            } else if (mat->tex_spec) {

                // reflective material
                v3 sample = img_sample(mat->tex_spec, rch.uv.x, rch.uv.y);

                // create a reflective ray
                ray_t ray_refl = (ray_t){ v3_add(rch.point, v3_scale(rch.normal, E2)), v3_refl(cur.ray.dir, rch.normal) };

                // add a render item
                rstk[++rstk_i] = (struct RenderItem){ cur.depth+1, ray_refl, v3_mul(cur.mask, sample) };
            } else if (mat->tex_refr) {
                // just refractive 


                // refractive material
                v3 sample = img_sample(mat->tex_refr, rch.uv.x, rch.uv.y);


                // create a refraction ray
                ray_t ray_refr = (ray_t){ 
                    v3_add(rch.point, v3_scale(rch.normal, rch.is_inside ? -E2 : E2)), 
                    v3_neg(v3_refract(cur.ray.dir, rch.normal, mat->IOR))
                };


                if (!v3_eqe(ray_refr.dir, V3_0, E2)) {
                //if (rch.is_inside) printf("is_inside: %d\n", (int)rch.is_inside);

                    // not total internal refraction
    //            printf("SDF: %f,%f,%f\n", V3__(ray_refr.dir));

                    // add a render item

                }
                rstk[++rstk_i] = (struct RenderItem){ cur.depth+1, ray_refr, V3_1 };
            }


        } else {
            // nothing?
        }

    }


/*

    rch_t cur = raycast(scene, ray);

    if (cur.hit && cur.hit->model.mat) {
        // we hit something, return its color

        mat_t mat = *cur.hit->model.mat;

        if (mat.tex_diff.data) {

            //return mat.tex_diff.data[0];
        } else {

            // for right now, just do this
            float val = v3_dot(cur.normal, v3_add(V3_DOWN, V3_FORWARD));
            val = val > 0 ? 0.0f : -val;

            // right now, just a red color
            return PIX_RGBA(val * 100.0f, 0, 0, 0);
        }


    } else {
        // nothing was hit, return a background color
        return PIX_RGBA(32, 32, 32, 255);

    }*/

    return vsum;

}

// main render function
void cdrt_render(img_t dest, scene_t scene) {
    // recalculate the cache
    cdrt_scene_recache(scene);

    obj_t cam = scene->cache.cam;

    if (cam == NULL) {
        fprintf(stderr, "CDRT: No camera!\n");
        return;
    }

    int i;
    // set image to black
    for (i = 0; i < dest->w * dest->h; ++i) dest->data[i] = PIX_RGBA(0, 0, 0, 0);

    // precompute this
    float tanHalfFOV = tanf(3.1415926535f * cam->camera.FOV * 0.5f / 180.0f);
    float aspectRatio = (float)dest->w / dest->h;

    // go through rows
    #pragma omp parallel for
    for (i = 0; i < dest->h; ++i) {
        int j;
        // go through columns
        for (j = 0; j < dest->w; ++j) {
            // compute local direction
            v3 local_dir = V3(
                tanHalfFOV * (2.0f * j / dest->w - 1.0f),
                (tanHalfFOV / aspectRatio) * (1.0f - 2.0f * i / dest->h),
                1.0f
            );

            // construct a ray from it, in global space
            ray_t view_ray = (ray_t){ 
                TRANSFORM_POS(cam->cache.baked_T), 
                //v3_unit(local_dir),
                v3_unit(m4T_transform_direction(cam->cache.baked_T, local_dir)) 
            };

            // calculate a pixel value
            v3 res = render_ray(scene, view_ray);

            // clamp values
            if (res.x < 0) res.x = 0; else if (res.x > 1) res.x = 1;
            if (res.y < 0) res.y = 0; else if (res.y > 1) res.y = 1;
            if (res.z < 0) res.z = 0; else if (res.z > 1) res.z = 1;

            // set the pixel data
            dest->data[i * dest->w + j] = PIX_RGBA(255.f * res.x, 255.f * res.y, 255.f * res.z, 255);
        }
    }

    //printf(M4X4_FMT "\n", M4X4__(cam->cache.baked_T));

    return;
}

