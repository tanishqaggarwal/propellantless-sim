#ifndef PSIM_CONTROLLER_TORQUE_HPP_
#define PSIM_CONTROLLER_TORQUE_HPP_

#include <psim/controller/torque.yml.hpp>

namespace psim {

class SatTorqueController : public SatTorqueControllerInterface<SatTorqueController> {
 private:
  typedef SatTorqueControllerInterface<SatTorqueController> Super;

 public:
  using Super::SatTorqueControllerInterface;

  SatTorqueController() = delete;
  virtual ~SatTorqueController() = default;

  /**
   * Lazy state fields.
   */
  Vector3 control_torque() const;
};
}  // namespace psim

#endif
