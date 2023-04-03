#pragma once
#include <functional>
#include "data.hpp"

class Weno : public Data {

protected:
  int k; //polynomial degree

  const float eps = 1e-6; //for alpha coeffs

  std::vector<std::vector<float>> c_rj; //polynomial coeffs for degree k-1
  std::vector<float> d;

  float ComputePolynomialCoeff(int r, int j);

  void ComputeBetas(std::vector<float>& dst, int i);

  void ComputeWeightedAlphas(
      std::vector<float>& dst, int i);

  void ComputeWenoApproximation(std::vector<float>& fluxes);

  void RungeKuttaTimeStep(
    std::function<float(float, float)> L);

  Weno(const nlohmann::json& cfg) : Data(cfg) {
    k = cfg["k"];

    ge_w = k;
    mesh.resize(int(cfg["nodes_num"]) + 2 * ge_w);

    //compute polynomial coeffs
    c_rj.resize(k + 1);
    for (int r = -1; r < k; r++) {
      c_rj[r + 1].resize(k);
      for (int j = 0; j < k; ++j) {
          c_rj[r + 1][j] = ComputePolynomialCoeff(r, j);
      }
    }

    for (int i = 0; i < c_rj.size(); ++i) {
      for (int j = 0; j < c_rj[i].size(); ++j) {
        std::cout << c_rj[i][j] << " ";
      }
      std::cout << std::endl;
    }

    //initialize d coeffs
    switch (k) {
    case 1:
      d.push_back(1.0);
      break;
    case 2:
      d.push_back(2.0 / 3.0);
      d.push_back(1.0 / 3.0);
      break;
    case 3:
      d.push_back(0.3);
      d.push_back(0.6);
      d.push_back(0.1);
      break;
    default:
      throw std::runtime_error (
        "Sorry, no d coeffs for k > 3 yet :( ");
      break;
    }
  }
};

float Weno::ComputePolynomialCoeff(int r, int j) {
  //TODO: r = -1 not needed?
  float res = 0.0;
  for (int m = j + 1; m < k + 1; ++m) {
    float sum_num = 0.0;
    for (int l = 0; l < k + 1; ++l) {
      if (l == m)
        continue;
      float prod_num = 1.0;
      for (int q = 0; q < k + 1; ++q) {
        if (q == l || q == m)
          continue;
        prod_num *= r - q + 1;
      }
      sum_num += prod_num;
    }
    float prod_denom = 1.0;
    for (int l = 0; l < k + 1; ++l) {
      if (l == m)
        continue;
      prod_denom *= (m - l);
    }
    res += (sum_num / prod_denom);
  }
  return res;
}

inline float sqr(float x) {
  return x*x;
}

void Weno::ComputeBetas(
    std::vector<float> &dst, int i) {

  switch (k) {
  case 2:
    dst[0] = (mesh[i+1] - mesh[i])*(mesh[i+1] - mesh[i]);
    dst[1] = (mesh[i] - mesh[i-1])*(mesh[i] - mesh[i-1]);
    break;
  case 3:
    dst[0] = 13.0/12.0*sqr(mesh[i] - 2*mesh[i+1] + mesh[i+2]) +
        0.25*sqr(3*mesh[i] - 4*mesh[i+1] + mesh[i+2]);
    dst[1] = 13.0/12.0*sqr(mesh[i-1] - 2*mesh[i] + mesh[i+1])+
        0.25*sqr(mesh[i+1] - mesh[i-1]);
    dst[2] = 13.0/12.0*sqr(mesh[i-2] - 2*mesh[i-1] + mesh[i]) +
        0.25*sqr(mesh[i-2] - 4*mesh[i-1] + 3*mesh[i]);
    break;
  default:
    throw std::runtime_error (
      "Sorry, beta coeffs k = 2 or 3 only :( ");
    break;
  }

  return;

}

void Weno::ComputeWeightedAlphas(
    std::vector<float> &dst, int i) {
  std::vector<float> betas(k);
  ComputeBetas(betas, i);
  for (int r = 0; r < k; ++r) {
    dst[r] = d[r] / (eps + betas[r]) / (eps + betas[r]);
  }
  float sum = 0.f;
  for (float a : dst) {
    sum += a;
  }
  for (float& a : dst) {
    a /= sum;
  }
  return;
}

void Weno::ComputeWenoApproximation(std::vector<float>& fluxes) {

  fluxes.clear();
  fluxes.resize(mesh.size());

  for (int i = ge_w; i < mesh.size() - ge_w; ++i) {
    std::vector<float> recs_in_stencils(k); //function reconstruction for different stencils
    for (int r = 0; r < k; ++r) { //different stencils for i-th node
      for (int j = 0; j < k; ++j) {
        recs_in_stencils[r] += c_rj[r + 1][j] * mesh[i - r + j];
      }
    }
    std::vector<float> w(k);
    ComputeWeightedAlphas(w, i);
    for (int r = 0; r < k; ++r) {//k stencils (biased to the left)
      fluxes[i] += w[r] * recs_in_stencils[r];
    }
  }

  return;
}

void Weno::RungeKuttaTimeStep(
  std::function<float(float, float)> L) {

  std::vector<float> fluxes;
  ComputeWenoApproximation(fluxes);

  std::vector<float> u_n(mesh);

  std::vector<float> u1(u_n.size());
  for (size_t i = ge_w; i < u_n.size() - ge_w - 1; ++i) {
    u1[i] = u_n[i] + delta_t * L(fluxes[i + 1], fluxes[i]);
  }
  mesh = u1;
  if (time_ord == 1) {
    return;
  }

  //mesh is now u1
  bc->ApplyBoundaryConditions(this);
  ComputeWenoApproximation(fluxes);
  std::vector<float> u2(u_n.size());
  for (size_t i = ge_w; i < u_n.size() - ge_w - 1; ++i) {
    u2[i] = 0.75*u_n[i] + 0.25 * u1[i] + 0.25 * delta_t * L(fluxes[i + 1], fluxes[i]);
  }
  mesh = u2;
  if (time_ord == 2) {
    return;
  }

  //mesh is now u2
  bc->ApplyBoundaryConditions(this);
  ComputeWenoApproximation(fluxes);
  for (size_t i = ge_w; i < u_n.size() - ge_w - 1; ++i) {
    mesh[i] = 0.3333*u_n[i] + 0.6667 * u2[i] + 0.6667 * delta_t * L(fluxes[i + 1], fluxes[i]);
  }

  return;
}
