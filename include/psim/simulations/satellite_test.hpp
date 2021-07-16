#ifndef PSIM_SIMULATIONS_SATELLITE_TEST_HPP_
#define PSIM_SIMULATIONS_SATELLITE_TEST_HPP_

#include <psim/core/configuration.hpp>
#include <psim/core/model_list.hpp>

namespace psim {

/** @brief Models a single satellite's attitude and orbital dynamics. All models
 *         are backed by flight software's GNC implementations if possible.
 */
class TwoPointSatelliteTest : public ModelList {
 public:
  TwoPointSatelliteTest() = delete;
  virtual ~TwoPointSatelliteTest() = default;

  TwoPointSatelliteTest(RandomsGenerator &randoms, Configuration const &config);
};
}  // namespace psim

#endif