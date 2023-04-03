#pragma once
#include "init.hpp"

class BoundaryConditionsApplier {
public:
  virtual void ApplyBoundaryConditions(Data* dst) = 0;
};


class Symmetry : public BoundaryConditionsApplier {
public:
  void ApplyBoundaryConditions(Data* dst) override;
};

class Periodic : public BoundaryConditionsApplier {
public:
  void ApplyBoundaryConditions(Data* dst) override;
};
