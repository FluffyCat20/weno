#pragma once
#include <vector>
#include <nlohmann/json.hpp>
#include <fstream>
#include <memory>

class Data;

class Initializer {

public:
  virtual void Init(Data* dst) = 0;
};

class SimpleDiscontiniuty : public Initializer {
  float v0_, v1_; //values on left and right from discontiniuty
  float x0_; //discontiniuty coordinate (from 0 to 1)

public:

  SimpleDiscontiniuty(float x0, float v0, float v1) :
    x0_(x0), v0_(v0), v1_(v1){};

  SimpleDiscontiniuty(const nlohmann::json& cfg) {
    x0_ = cfg["x0"];
    v0_ = cfg["v0"];
    v1_ = cfg["v1"];
  }

  void Init(Data* dst) override;
};

class JiangShuFirstTest : public Initializer {
  float a = 0.5;
  float z = -0.7;
  float delta = 0.005;
  float alpha = 10.0;
  float beta = log(2.0) / 36 / delta / delta;

  float G(float x, float y);
  float F(float x, float y);

public:
  JiangShuFirstTest(){}
  void Init(Data* dst) override;
};

float JiangShuFirstTest::G(float x, float y) {
  float pow = -beta * (x - y) * (x - y);
  return exp(pow);
}

float JiangShuFirstTest::F(float x, float y) {
  return std::sqrt(std::max(0.0,
         1.0 - alpha*alpha*(x - y)*(x - y)));
}
