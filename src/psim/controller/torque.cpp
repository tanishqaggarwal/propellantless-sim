#include <psim/controller/torque.hpp>

#include <gnc/utilities.hpp>

#include <lin/core.hpp>
#include <lin/generators/constants.hpp>
#include <lin/math.hpp>
#include <lin/references.hpp>

namespace psim {

Vector3 SatTorqueController::control_torque() const
{
  return Vector3();
}

}