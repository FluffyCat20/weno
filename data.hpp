#pragma once
#include <iostream>
#include <stdexcept>
#include "boundary_conditions.hpp"

class OutputManager;

class Data {

protected:
  std::vector<float> mesh;
  float l, r; //area bounds coordinates
  float delta_x;
  int ge_w; //ghost edges width
  float current_t = 0.f;
  float delta_t = 0.f;
  float courant_number = 0.0;
  float t_end = 0.0;
  int time_ord = 0;

  friend class SimpleDiscontiniuty;
  friend class JiangShuFirstTest;

  friend class Symmetry;
  friend class Periodic;

  friend class OutputManager;

  virtual void CalcDeltaT() = 0;

public:

  Data(const nlohmann::json& cfg) {
    r = float(cfg["right_bound"]);
    l = float(cfg["left_bound"]);
    delta_x = (r - l) /
        float(cfg["nodes_num"]);
    current_t = 0.f;
    courant_number = float(cfg["courant_num"]);
    t_end = float(cfg["t_end"]);
    time_ord = float(cfg["time_order"]);
  }

  virtual void DoStep() = 0;

  std::shared_ptr<BoundaryConditionsApplier> bc;

  void print() {
    for (auto val : mesh) {
      std::cout << val << " " << " ";
    }
    std::cout << std::endl;
  }

  const float CurrentTime() {
    return current_t;
  }

  const float TimeEnd() {
    return t_end;
  }

};

void SimpleDiscontiniuty::Init(Data* dst) {
  size_t disc_node = size_t(dst->mesh.size() * x0_);
  for (size_t i = dst->ge_w; i < disc_node; ++i) {
    dst->mesh[i] = v0_;
  }
  for (size_t i = disc_node; i < dst->mesh.size() - dst->ge_w; ++i) {
    dst->mesh[i] = v1_;
  }
  return;
}

void JiangShuFirstTest::Init(Data *dst) {
  for (int i = dst->ge_w; i < dst->mesh.size() - dst->ge_w; ++i) {
    float x = dst->l + i * dst->delta_x;
    if (-0.8 <= x && x <= -0.6) {
      float G1 = G(x, z - delta);
      float G2 = G(x, z + delta);
      float G3 = G(x, z);
      dst->mesh[i] = (G1 + G2 + 4*G3)/6.0;
    }
    else if (-0.4 <= x && x <= -0.2) {
      dst->mesh[i] = 1.0;
    }
    else if (0.0 <= x && x <= 0.2) {
      dst->mesh[i] = 1.0 - 10.0 * std::abs(x - 0.1);
    }
    else if (0.4 <= x && x <= 0.6) {
      float F1 = F(x, a - delta);
      float F2 = F(x , a + delta);
      float F3 = F(x, a);
      dst->mesh[i] = (F1 + F2 + 4*F3)/6.0;
    }
    else {
      dst->mesh[i] = 0.0;
    }
  }
}

void Symmetry::ApplyBoundaryConditions(Data* dst) {
  for (size_t i = 1; i < dst->ge_w + 1; ++i) {
    dst->mesh[dst->ge_w - i] = dst->mesh[dst->ge_w + i - 1];
    dst->mesh[dst->mesh.size() - dst->ge_w - 1 + i] =
        dst->mesh[dst->mesh.size() - dst->ge_w - i];
  }
  return;
}

void Periodic::ApplyBoundaryConditions(Data *dst) {
  for (size_t i = 0; i < dst->ge_w; ++i) {
    size_t right_pos = dst->mesh.size() - 2*dst->ge_w + i;
    dst->mesh[i] = dst->mesh[right_pos];
    //std::cout << right_pos << " copied to " << i << std::endl;
    dst->mesh[dst->mesh.size() - dst->ge_w + i] = dst->mesh[dst->ge_w + i];
    //std::cout << dst->ge_w + i << " copied to "
      //<< dst->mesh.size() - dst->ge_w + i << std::endl;
  }
  return;
}
