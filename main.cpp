#include <cstring>
#include "WenoSimpleTransfer.hpp"
#include "printer.hpp"

int main(int argc, char* argv[])
{
  if (argc < 3 || strcmp(argv[1], "-c")) {
    throw std::runtime_error(
      "Usage: -c config_path");
  }

  std::ifstream input(argv[2]);
  if (input.fail()) {
    throw std::runtime_error("config not found");
  }

  nlohmann::json config_data;
  input >> config_data;

  std::cout << config_data << std::endl;

  std::shared_ptr<Data> resolver;
  if (config_data["type"] == "WenoSimpleTransfer") {
    resolver.reset(new WenoSimpleTransfer(config_data));
  } else {
    std::cout << "Unkown resolver type" << std::endl;
    return 0;
  }

  std::shared_ptr<Initializer> initializer;
  std::string init_type = config_data["init_type"];
  if (init_type == "simple_discontinuity") {
    initializer.reset(new SimpleDiscontiniuty(config_data));
  }
  else if (init_type == "jiang_shu_first_test") {
    initializer.reset(new JiangShuFirstTest());
  } else {
    throw std::runtime_error("unknown initial conditions type :(");
  }

  auto& bc_applier = resolver->bc;
  std::string bc_type = config_data["boundary_conditions_type"];
  if (bc_type == "symmetry") {
    bc_applier.reset(new Symmetry());
  }
  else if (bc_type == "periodic") {
    bc_applier.reset(new Periodic());
  } else {
    throw std::runtime_error("unknown boundary conditions type :(");
  }

  initializer->Init(resolver.get());
  bc_applier->ApplyBoundaryConditions(resolver.get());

  OutputManager printer(config_data, resolver.get());

  while (resolver->CurrentTime() < resolver->TimeEnd()) {
    resolver->DoStep();
    bc_applier->ApplyBoundaryConditions(resolver.get());
    printer.SaveMeshToFile(resolver.get());
  }

  return 0;
}
