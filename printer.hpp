#pragma once
#include <iomanip>
#include "data.hpp"

class OutputManager {
  std::ofstream mesh_outfile;

public:

  OutputManager(const nlohmann::json& cfg, const Data* src) {
    mesh_outfile.open(cfg["outfile"]);
    if (!mesh_outfile.good()) {
      throw (std::runtime_error("Cannot open output file"));
    }
    mesh_outfile << "VARIABLES = \"x\", \"v\"" << std::endl;
    mesh_outfile << std::setprecision(4);
    SaveMeshToFile(src);
  }

  ~OutputManager() {
    mesh_outfile.close();
  }

  void SaveMeshToFile(const Data* src);
};

void OutputManager::SaveMeshToFile(const Data* src) {
  mesh_outfile <<
    "ZONE T=\"timezone\" SOLUTIONTIME=" <<
    src->current_t << std::endl <<
    "I = " << src->mesh.size() - src->ge_w * 2 <<
    ", F=POINT" << std::endl;
  for (size_t i = src->ge_w; i < src->mesh.size() - src->ge_w; ++i) {
    mesh_outfile << src->l + (i - src->ge_w + 1) * src->delta_x << " " <<
                    src->mesh[i] << std::endl;
  }

  std::cout << "Done: time = " << src->current_t << std::endl;

  return;
}
