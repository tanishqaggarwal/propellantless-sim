#include <psim/simulations/satellite_test.hpp>
#include <psim/controller/torque.hpp>
#include <psim/truth/earth.hpp>
#include <psim/truth/satellite.hpp>
#include <psim/truth/time.hpp>

namespace psim {

TwoPointSatelliteTest::TwoPointSatelliteTest(
    RandomsGenerator &randoms, Configuration const &config)
  : ModelList(randoms) {
  add<Time>(randoms, config);
  add<EarthGnc>(randoms, config);
  add<SatTorqueController>(randoms, config);
  add<TwoPointSatellite>(randoms, config);
}
} // namespace psim
