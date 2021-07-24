#include <psim/truth/satellite.hpp>

#include <gnc/utilities.hpp>

#include <lin/core.hpp>
#include <lin/generators/constants.hpp>
#include <lin/math.hpp>
#include <lin/references.hpp>

#include <psim/truth/attitude_utilities.hpp>
#include <psim/truth/orbit_utilities.hpp>

namespace psim {

static constexpr Real half = 0.5;

Vector3 TwoPointSatellite::qdot_to_body_w(Vector4 const& qdot_body_eci, Vector4 const& q_body_eci)
{
  Vector4 half_w, qc;
  gnc::utl::quat_conj(q_body_eci, qc);
  gnc::utl::quat_cross_mult(qc, qdot_body_eci, half_w);
  return 2*lin::ref<Vector3>(half_w, 0, 0);
}

Vector4 TwoPointSatellite::body_w_to_qdot(Vector3 const& w_body, Vector4 const& q_body_eci)
{
  Vector4 half_w, ret;
  lin::ref<Vector3>(half_w, 0, 0) = half*w_body;
  gnc::utl::quat_cross_mult(q_body_eci, half_w, ret);
  return ret;
}

Vector<TwoPointSatellite::N_ODE>
TwoPointSatellite::dynamics(Real t, Vector<TwoPointSatellite::N_ODE> const &x, void *ptr)
{
  /*
   * Step 0. Unpack and prepare values that are useful throughout the computation.
   */
  auto const *data = static_cast<TwoPointSatellite::IntegratorData *>(ptr);

  auto const &m_A = data->mA;
  auto const &m_B = data->mB;
  auto const &d = data->d;
  auto const &q_eci_ecef = data->q_eci_ecef;

  auto const r_eci = lin::ref<Vector3>(x, 0, 0);
  auto const v_eci = lin::ref<Vector3>(x, 3, 0);
  Vector4 q_body_eci = lin::ref<Vector4>(x, 6, 0);
  auto const w_body = lin::ref<Vector3>(x, 10, 0);

  Vector<N_ODE> dx;

  Real m = m_A + m_B;
  Vector4 q_ecef_eci, q_eci_body;
  gnc::utl::quat_conj(q_eci_ecef, q_ecef_eci);
  gnc::utl::quat_conj(q_body_eci, q_eci_body);

  /*
   * Step 1. Compute vdot.
   */

  /*
   * Step 1.1. Get locations of satellite components in ECEF.
   */
  Vector3 xA, xB;
  get_mass_locations_eci(r_eci, q_eci_body, d, m_A, m_B, xA, xB);
  gnc::utl::rotate_frame(q_ecef_eci, xA);
  gnc::utl::rotate_frame(q_ecef_eci, xB);

  /*
   * Step 1.2. Call GEOGRAV.
   */
  Vector3 F_a = m_A * orbit::gravity(xA);
  Vector3 F_b = m_B * orbit::gravity(xB);

  /*
   * Step 1.3. Transform gravity vector into ECI.
   */
  gnc::utl::rotate_frame(q_eci_ecef, F_a);
  gnc::utl::rotate_frame(q_eci_ecef, F_b);

  /*
   * Step 1.4. Sum forces and compute total force.
   */
  Vector3 a_eci = (F_a + F_b) / m;

  /*
   * Step 2. Compute wdot.
   */
  
  /*
   * Step 2.1. Transform gravity forces into body frame.
   */
  gnc::utl::rotate_frame(q_body_eci, F_a);
  gnc::utl::rotate_frame(q_body_eci, F_b);

  /*
   * Step 2.2. Compute rotational inertia.
   */
  get_mass_locations_body(d, m_A, m_B, xA, xB);
  Real rA = lin::norm(xA), rB = lin::norm(xB);
  Real J = rA*rA*m_A + rB*rB*m_B;

  /*
   * Step 2.3. Compute angular acceleration.
   */
  Vector3 alpha_body = (lin::cross(xA, F_a) + lin::cross(xB, F_b)) / J;

  /*
   * Step 3. Compute change in quaternion.
   */
  Vector4 qdot_body_eci = body_w_to_qdot(w_body, q_body_eci);

  /**
   * Step 4. Populate the derivative vector
   */
  lin::ref<Vector3>(dx, 0, 0) = v_eci;
  lin::ref<Vector3>(dx, 3, 0) = a_eci;
  lin::ref<Vector4>(dx, 6, 0) = qdot_body_eci;
  lin::ref<Vector3>(dx, 10, 0) = alpha_body;

  return dx;
}

void TwoPointSatellite::get_mass_locations_body(Real l, Real mA, Real mB,
  Vector3& xA_body, Vector3& xB_body)
{
  Real m = mA + mB;
  Vector3 d; d(0) = l;
  xA_body = - d * mB / m;
  xB_body = d * mA / m;
}

void TwoPointSatellite::get_mass_locations_eci(Vector3 const& x_eci,
  Vector4 const& q_eci_body, Real l, Real mA, Real mB, Vector3& xA_eci,
  Vector3& xB_eci)
{
  get_mass_locations_body(l, mA, mB, xA_eci, xB_eci);
  
  gnc::utl::rotate_frame(q_eci_body, xA_eci);
  gnc::utl::rotate_frame(q_eci_body, xB_eci);
  xA_eci = xA_eci + x_eci;
  xB_eci = xB_eci + x_eci;
}

void TwoPointSatellite::get_mass_locations_eci(Vector3& xA_eci,
  Vector3& xB_eci) const
{
  auto const &x_eci = truth_satellite_orbit_r.get();
  auto const &l = truth_satellite_length.get();
  auto const &q_body_eci = truth_satellite_attitude_q_body_eci.get();
  Real mA = truth_satellite_A_m.get();
  Real mB = truth_satellite_B_m.get();

  Vector4 q_eci_body;
  gnc::utl::quat_conj(q_body_eci, q_eci_body);
  get_mass_locations_eci(x_eci, q_eci_body, l, mA, mB, xA_eci, xB_eci);
}

void TwoPointSatellite::step()
{
  this->Super::step();

  auto const &dt = truth_dt_s->get();
  auto const &q_eci_ecef = truth_earth_q_eci_ecef->get();
  auto const &mA = truth_satellite_A_m.get();
  auto const &mB = truth_satellite_B_m.get();
  auto const &d = truth_satellite_length.get();

  auto &r_eci = truth_satellite_orbit_r.get();
  auto &v_eci = truth_satellite_orbit_v.get();
  auto &q_body_eci = truth_satellite_attitude_q_body_eci.get();
  auto &w_body = truth_satellite_attitude_w.get();

  IntegratorData data{mA, mB, d, q_eci_ecef};

  // Prepare integrator inputs.
  Vector<N_ODE> x;
  lin::ref<Vector3>(x, 0, 0) = r_eci;
  lin::ref<Vector3>(x, 3, 0) = v_eci;
  lin::ref<Vector4>(x, 6, 0) = q_body_eci;
  lin::ref<Vector3>(x, 10, 0) = w_body;

  // Simulate dynamics.
  x = ode(Real(0.0), dt, x, &data, TwoPointSatellite::dynamics);

  // Write back to our state fields
  r_eci = lin::ref<Vector3>(x, 0, 0);
  v_eci = lin::ref<Vector3>(x, 3, 0);
  q_body_eci = lin::ref<Vector4>(x, 6, 0);
  w_body = lin::ref<Vector3>(x, 10, 0);
}

Vector3 TwoPointSatellite::truth_satellite_orbit_A_r() const
{
  Vector3 xA, xB;
  get_mass_locations_eci(xA, xB);
  return xA;
}

Vector3 TwoPointSatellite::truth_satellite_orbit_B_r() const
{
  Vector3 xA, xB;
  get_mass_locations_eci(xA, xB);
  return xB;
}

Real TwoPointSatellite::truth_satellite_orbit_T() const
{
  auto const &v = truth_satellite_orbit_v.get();
  auto const &w = truth_satellite_attitude_w.get();
  Real mA = truth_satellite_A_m.get();
  Real mB = truth_satellite_B_m.get();

  // Compute twice the translational kinetic energy
  Real m = mA + mB;
  Real s = lin::norm(v);
  Real tKE = m * s * s;

  // Compute rotational kinetic energy by first:
  // 1. Computing distances of masses from COM
  // 2. Computing moment of inertia J
  // 3. Computing J w^2
  Vector3 xA, xB;
  get_mass_locations_eci(xA, xB);
  Real rA = lin::norm(xA);
  Real rB = lin::norm(xB);
  Real J = rA*rA*mA + rB*rB*mB;
  Real ws = lin::norm(w);
  Real rKE = J * ws * ws;

  return half * (tKE + rKE);
}

Real TwoPointSatellite::truth_satellite_orbit_U() const
{
  Vector3 xA, xB;
  get_mass_locations_eci(xA, xB);
  Real rA = lin::norm(xA);
  Real rB = lin::norm(xB);
  Real mA = truth_satellite_A_m.get();
  Real mB = truth_satellite_B_m.get();

  return gnc::constant::mu_earth * (mA / rA + mB / rB);
}

Real TwoPointSatellite::truth_satellite_orbit_E() const
{
  Real T = truth_satellite_orbit_T();
  Real U = truth_satellite_orbit_U();

  return T - U;
}

}

