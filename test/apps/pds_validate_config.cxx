#include "conffwk/Configuration.hpp"
#include "confmodel/Session.hpp"
#include "confmodel/DaqApplication.hpp"
#include "appmodel/ConfigurationHelper.hpp"
#include "appmodel/SmartDaqApplication.hpp"
#include <iostream>
#include <memory>

int main(int argc, char** argv)
{
  if (argc < 4) return 2;
  try {
    for (int i = 3; i < argc; ++i) {
      // Module generation modifies the database; use a fresh view for each app.
      dunedaq::conffwk::Configuration db(std::string("oksconflibs:") + argv[1]);
      const auto* session = db.get<dunedaq::confmodel::Session>(argv[2]);
      if (!session) return 3;
      auto helper = std::make_shared<dunedaq::appmodel::ConfigurationHelper>(session);
      const auto* app = db.get<dunedaq::confmodel::DaqApplication>(argv[i]);
      if (!app) return 4;
      if (const auto* smart = app->cast<dunedaq::appmodel::SmartDaqApplication>()) {
        smart->generate_modules(helper);
      }
      std::cout << "VALIDATED " << argv[i] << " " << app->get_modules().size() << " modules\n";
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
