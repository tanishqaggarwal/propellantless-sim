#ifndef PSIM_TRUTH_SATELLITE_HPP_
#define PSIM_TRUTH_SATELLITE_HPP_

#include <psim/truth/satellite.yml.hpp>
#include <gnc/ode4.hpp>

namespace psim {

class TwoPointSatellite : public TwoPointSatelliteInterface<TwoPointSatellite> {
 private:
  typedef TwoPointSatelliteInterface<TwoPointSatellite> Super;

 protected:
  /**
   * Gets the locations of the two satellite components in body frame.
   * 
   * @param[in] l Length of satellite.
   * @param[in] mA Mass of satellite component A.
   * @param[in] mB Mass of satellite component B.
   * @param[out] xA_body Location of component A, in body frame.
   * @param[out] xB_body Location of component B, in body frame.
   */
  static void get_mass_locations_body(Real l, Real mA, Real mB,
    Vector3& xA_body, Vector3& xB_body);

  /**
   * Gets the locations of the two satellite components in ECI frame.
   * 
   * @param[in] x_eci Location of satellite COM in ECI.
   * @param[in] q_eci_body Quaternion transforming vectors from body frame to ECI.
   * @param[in] l Length of satellite.
   * @param[in] mA Mass of satellite component A.
   * @param[in] mB Mass of satellite component B.
   * @param[out] xA_eci Location of component A, in ECI.
   * @param[out] xB_eci Location of component B, in ECI.
   */
  static void get_mass_locations_eci(Vector3 const& x_eci,
    Vector4 const& q_eci_body, Real l, Real mA, Real mB, Vector3& xA_eci,
    Vector3& xB_eci);

  /**
   * Gets the locations of the two satellite components in ECI frame.
   * This function uses the state-field values of r, d, m_A, m_B and q in order
   * to compute this.
   * 
   * @param[out] xA Location of component A, in ECI.
   * @param[out] xB Location of component B, in ECI.
   */
  void get_mass_locations_eci(Vector3& xA_eci, Vector3& xB_eci) const;

  /**
   * Utility function for conversion between vector angular rate
   * and a quaternion representing the rotation rate.
   * 
   * @param w_body Angular rate in body frame
   * @param q_body_eci Quaternion rotating from ECI to body
   * @return Derivative of quaternion in ECI frame.
   */
  static Vector4 body_w_to_qdot(Vector3 const& w_body,
    Vector4 const& q_body_eci);

  /**
   * Utility function for conversion between vector angular rate
   * and a quaternion representing the rotation rate.
   * 
   * @param qdot_body_eci Derivative of quaternion that rotates from ECI to
   *                      body frame
   * @param q_body_eci Quaternion rotating from ECI to body frame
   * @return Angular rate vector in body frame
   */
  static Vector3 qdot_to_body_w(Vector4 const& qdot_body_eci,
    Vector4 const& q_body_eci);

  struct IntegratorData
  {
    Real const &mA;
    Real const &mB;
    Real const &d;
    Vector4 const &q_eci_ecef;
  };

  // Number of components in dynamics solver
  static constexpr size_t N_ODE = 13;

  // ODE integrator
  gnc::Ode4<Real, N_ODE> ode;

  /**
   * Applies full dual-point mass dynamics to the spacecraft.
   */
  static Vector<N_ODE> dynamics(Real t, Vector<N_ODE> const &x, void *ptr);

 public:
  using Super::TwoPointSatelliteInterface;

  TwoPointSatellite() = delete;
  virtual ~TwoPointSatellite() = default;

  virtual void step() override;

  /**
   * Lazy state fields.
   */

  Vector3 truth_satellite_orbit_A_r() const;
  Vector3 truth_satellite_orbit_B_r() const;
  Real truth_satellite_orbit_T() const;
  Real truth_satellite_orbit_U() const;
  Real truth_satellite_orbit_E() const;
};
}  // namespace psim

#endif