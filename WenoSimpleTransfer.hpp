#pragma once
#include "weno.hpp"

class WenoSimpleTransfer : public Weno {
  float a; //coef from equation

  void CalcDeltaT() override;

public:

  WenoSimpleTransfer(const nlohmann::json& cfg)
    : Weno(cfg) {
    a = cfg["a"];
  }

  void DoStep() override;

};

void WenoSimpleTransfer::CalcDeltaT() {
  delta_t = delta_x * courant_number / std::abs(a);
  if (current_t + delta_t > t_end) {
    delta_t = t_end - current_t;
  }
  current_t += delta_t;
  return;
}

void WenoSimpleTransfer::DoStep() {
  CalcDeltaT();
  RungeKuttaTimeStep([this](float f1, float f2) -> float {
    return -a / delta_x * (f1 - f2);
  });

  return;
}
