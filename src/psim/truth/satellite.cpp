#include <psim/truth/satellite.hpp>

#include <gnc/utilities.hpp>

#include <lin/core.hpp>
#include <lin/generators/constants.hpp>
#include <lin/math.hpp>
#include <lin/references.hpp>

#include <psim/truth/attitude_utilities.hpp>
#include <psim/truth/orbit_utilities.hpp>

namespace psim {

void TwoPointSatellite::get_mass_locations(Vector3& xA, Vector3& xB) const {
  auto const &x = truth_satellite_orbit_r.get();
  auto const &l = truth_satellite_length.get();
  auto const &q_body_eci = truth_satellite_attitude_q_body_eci.get();
  Real mA = truth_satellite_A_m.get();
  Real mB = truth_satellite_B_m.get();
  Real m = mA + mB;
  
  Vector3 d; // Vector from mass A to mass B
  d[0] = l;
  gnc::utl::rotate_frame(q_body_eci, d);

  xA = x - d * mB / m;
  xB = x + d * mA / m;
}

void TwoPointSatellite::step() {
  this->Super::step();

  struct IntegratorData
  {
    Real const &m;
    Vector3 const &earth_w;
    Vector3 const &earth_w_dot;
  };

  auto const &dt = truth_dt_s->get();
  auto const &earth_w = truth_earth_w->get();
  auto const &earth_w_dot = truth_earth_w_dot->get();
  auto const &q_eci_ecef = truth_earth_q_eci_ecef->get();
  auto const &mA = truth_satellite_A_m.get();
  auto const &mB = truth_satellite_B_m.get();
  Real m = mA + mB;

  auto &r_eci = truth_satellite_orbit_r.get();
  auto &v_eci = truth_satellite_orbit_v.get();
  auto &q_body_eci = truth_satellite_attitude_q_body_eci.get();
  auto &w_body = truth_satellite_attitude_w.get();

  // Rotate vectors into ECEF: required for geograv
  Vector4 q_ecef_eci;
  Vector3 r_ecef, v_ecef;
  gnc::utl::quat_conj(q_eci_ecef, q_ecef_eci);
  gnc::utl::rotate_frame(q_ecef_eci, r_eci, r_ecef);
  gnc::utl::rotate_frame(q_ecef_eci, v_eci, v_ecef);

  // Prepare integrator inputs.
  Vector<N_ODE> x;
  lin::ref<Vector3>(x, 0, 0) = r_ecef;
  lin::ref<Vector3>(x, 3, 0) = v_ecef;
  lin::ref<Vector4>(x, 6, 0) = q_body_eci;
  lin::ref<Vector3>(x, 10, 0) = w_body;
  IntegratorData data{m, earth_w, earth_w_dot};

  // Simulate dynamics.
  x = ode(Real(0.0), dt, x, &data,
      [](Real t, Vector<N_ODE> const &x, void *ptr) -> Vector<N_ODE> {
        auto const *data = static_cast<IntegratorData *>(ptr);

        auto const &m = data->m;
        auto const earth_w = (data->earth_w + t * data->earth_w_dot).eval();
        auto const &earth_w_dot = data->earth_w_dot;

        auto const r_ecef = lin::ref<Vector3>(x, 0, 0);
        auto const v_ecef = lin::ref<Vector3>(x, 3, 0);
        auto const q_body_eci = lin::ref<Vector4>(x, 6, 0);
        auto const w_body = lin::ref<Vector3>(x, 10, 0);

        Vector<N_ODE> dx;

        // Orbital dynamics
        {
          Vector3 const a_ecef = orbit::acceleration(
              earth_w, earth_w_dot, r_ecef.eval(), v_ecef.eval(), 0, m);

          lin::ref<Vector3>(dx, 0, 0) = v_ecef;
          lin::ref<Vector3>(dx, 3, 0) = a_ecef;
        }

        // Attitude dynamics - quaternion
        {
          Vector4 dq_body_eci;
          Vector4 const dq = {
              0.5 * w_body(0), 0.5 * w_body(1), 0.5 * w_body(2), 0.0};
          gnc::utl::quat_cross_mult(dq, q_body_eci.eval(), dq_body_eci);

          lin::ref<Vector4>(dx, 6, 0) = dq_body_eci;
        }

        return dx;
      });

  // Write back to our state fields
  r_ecef = lin::ref<Vector3>(x, 0, 0);
  v_ecef = lin::ref<Vector3>(x, 3, 0);
  gnc::utl::rotate_frame(q_eci_ecef, r_ecef, r_eci);
  gnc::utl::rotate_frame(q_eci_ecef, v_ecef, v_eci);

  q_body_eci = lin::ref<Vector4>(x, 6, 0);
  w_body = lin::ref<Vector3>(x, 10, 0);
}

Real TwoPointSatellite::truth_satellite_orbit_T() const {
  static constexpr Real half = 0.5;

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
  get_mass_locations(xA, xB);
  Real rA = lin::norm(xA);
  Real rB = lin::norm(xB);
  Real J = rA*rA*mA + rB*rB*mB;
  Real ws = lin::norm(w);
  Real rKE = J * ws * ws;

  return half * (tKE + rKE);
}

Real TwoPointSatellite::truth_satellite_orbit_U() const {
  Vector3 xA, xB;
  get_mass_locations(xA, xB);
  Real rA = lin::norm(xA);
  Real rB = lin::norm(xB);
  Real mA = truth_satellite_A_m.get();
  Real mB = truth_satellite_B_m.get();

  return gnc::constant::mu_earth * (mA / rA + mB / rB);
}

Real TwoPointSatellite::truth_satellite_orbit_E() const {
  Real T = truth_satellite_orbit_T();
  Real U = truth_satellite_orbit_U();

  return T - U;
}

}

