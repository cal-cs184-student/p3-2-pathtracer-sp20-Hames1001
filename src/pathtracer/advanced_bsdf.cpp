#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>


using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Implement MirrorBSDF
    *pdf = 1.0;
    reflect(wo, wi);
  return reflectance / abs_cos_theta(*wi);
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
  // TODO: proj3-2, part 3
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
    auto tan = sin_theta(h) / cos_theta(h);
    auto numerator = exp(-1 * (pow(tan, 2) / pow(alpha, 2)));
    auto denominator = PI * pow(alpha, 2) * pow(cos_theta(h), 4);
  return numerator / denominator;
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
  // TODO: proj3-2, part 3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Spectrum.
    
    auto eta_2k = eta * eta + k * k;
    auto nCos = 2.0 * eta * cos_theta(wi);
    auto R_s = (eta_2k - nCos + pow(cos_theta(wi), 2)) / (eta_2k + nCos + pow(cos_theta(wi), 2));
    auto R_p = ((eta_2k * pow(cos_theta(wi), 2)) - nCos + 1.0) / ((eta_2k * pow(cos_theta(wi), 2)) + nCos + 1.0);
  return (R_s + R_p) / 2.0;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO: proj3-2, part 3
  // Implement microfacet model here.
    Vector3D n = Vector3D(0,0,1);
    if (dot(n, wo) > 0 && dot(n, wi) > 0) {
        Vector3D h = wo + wi;
        h.normalize();
        return (F(wi) * G(wo, wi) * D(h)) / (4.0 * dot(n, wo) * dot(n, wi));
    }
    return Spectrum();
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: proj3-2, part 3
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.
  Vector2D sd = sampler.get_sample();
    auto rg1 = sd.x;
    auto rg2 = sd.y;
    auto thetaH = atan(sqrt(-1 * pow(alpha, 2) * log(1.0 - rg1)));
    auto phiH = 2.0 * PI * rg2;
    auto h = Vector3D(sin(thetaH) * sin(phiH), sin(thetaH) * cos(phiH), cos(thetaH));
    *wi = (2 * h * dot(h, wo)) - wo;
    //*wi = cosineHemisphereSampler.get_sample(pdf);
    if (wi->z <= 0) {
        *pdf = 0;
        return Spectrum();
    }
    auto p_theta = (2 * sin(thetaH) * exp(-1 * pow(tan(thetaH), 2) / pow(alpha, 2))) / (pow(alpha, 2) * pow(cos(thetaH), 3)) ;
    auto p_phi  = 1 / (2 * PI);
    auto pH = (p_theta * p_phi) / sin(thetaH);
    *pdf = pH / (4 * dot(*wi, h));
    return f(wo, *wi);
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Implement RefractionBSDF
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305
    bool inOrEx = refract(wo, wi, ior);
    if (!inOrEx) {
        reflect(wo, wi);
        *pdf = 1.0;
        return reflectance / abs_cos_theta(*wi);
    }
    auto R = (1.0 - ior) / (1.0 + ior);
    auto eta = ior;
    if (wo.z > 0) {
        eta = 1.0 / ior;
    }
    R = pow(R, 2);
    R = R + (1 - R) * pow(1 - abs_cos_theta(wo), 5);
    if (coin_flip(R)) {
        reflect(wo, wi);
        *pdf = R;
        return R * (reflectance / abs_cos_theta(*wi));
    } else {
        refract(wo, wi, ior);
        *pdf = 1 - R;
        return *pdf * (transmittance / abs_cos_theta(*wi) / pow(eta, 2));
    }
    
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
  // TODO:
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
    *wi = Vector3D(-wo.x, -wo.y, wo.z);
    
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {
  // TODO:
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
    int coordinateSign = wo.z > 0 ? -1 : 1;
    double n = ior;
    if (wo.z > 0) {
        n = 1.0 / ior;
    }
    auto sLaw = 1 - (pow(n, 2) * sin_theta2(wo));
    if (sLaw < 0) {
        return false;
    }
    *wi = Vector3D(-1 * n * wo.x, -1 * n * wo.y, coordinateSign * sqrt(sLaw));
  return true;
}

} // namespace CGL
