#ifndef PSIM_TRUTH_SATELLITE_HPP_
#define PSIM_TRUTH_SATELLITE_HPP_

#include <psim/truth/satellite.yml.hpp>
#include <gnc/ode4.hpp>

namespace psim {

class TwoPointSatellite : public TwoPointSatelliteInterface<TwoPointSatellite> {
 private:
  typedef TwoPointSatelliteInterface<TwoPointSatellite> Super;

 protected:
  static constexpr size_t N_ODE = 13; // Number of components in ODE
  gnc::Ode4<Real, N_ODE> ode;
  void get_mass_locations(Vector3& xA, Vector3& xB) const;

 public:
  using Super::TwoPointSatelliteInterface;

  TwoPointSatellite() = delete;
  virtual ~TwoPointSatellite() = default;

  virtual void step() override;

  Real truth_satellite_orbit_T() const;
  Real truth_satellite_orbit_U() const;
  Real truth_satellite_orbit_E() const;
};
}  // namespace psim

#endif