/*
MIT License

Copyright (c) 2019 Nathan Zimmerberg

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*
The MIT License (MIT); this license applies to GeographicLib,
versions 1.12 and later.

Copyright (c) 2008-2019, Charles Karney

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/

// ../geograv.hpp Generated by python script gravcodeupdate.py
/**
 * \file geograv.hpp
 * \author Nathan Zimmerberg
 * \date
 * \brief copied from SphericalEngine.cpp/hpp in https://geographiclib.sourceforge.io/ and modified to use static memory
 * The below mathy documentation is from https://geographiclib.sourceforge.io/ SphericalEngine.cpp:
 *
 * Copyright (c) Charles Karney (2011-2018) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 *
 * The general sum is\verbatim
 V(r, theta, lambda) = sum(n = 0..N) sum(m = 0..n)
   q^(n+1) * (C[n,m] * cos(m*lambda) + S[n,m] * sin(m*lambda)) * P[n,m](t)
\endverbatim
 * where <tt>t = cos(theta)</tt>, <tt>q = a/r</tt>.  In addition write <tt>u =
 * sin(theta)</tt>.
 *
 * <tt>P[n,m]</tt> is a normalized associated Legendre function of degree
 * <tt>n</tt> and order <tt>m</tt>.  Here the formulas are given for full
 * normalized functions (usually denoted <tt>Pbar</tt>).
 *
 * Rewrite outer sum\verbatim
 V(r, theta, lambda) = sum(m = 0..N) * P[m,m](t) * q^(m+1) *
    [Sc[m] * cos(m*lambda) + Ss[m] * sin(m*lambda)]
\endverbatim
 * where the inner sums are\verbatim
   Sc[m] = sum(n = m..N) q^(n-m) * C[n,m] * P[n,m](t)/P[m,m](t)
   Ss[m] = sum(n = m..N) q^(n-m) * S[n,m] * P[n,m](t)/P[m,m](t)
\endverbatim
 * Evaluate sums via Clenshaw method.  The overall framework is similar to
 * Deakin with the following changes:
 * - Clenshaw summation is used to roll the computation of
 *   <tt>cos(m*lambda)</tt> and <tt>sin(m*lambda)</tt> into the evaluation of
 *   the outer sum (rather than independently computing an array of these
 *   trigonometric terms).
 * - Scale the coefficients to guard against overflow when <tt>N</tt> is large.
 * .
 * For the general framework of Clenshaw, see
 * http://mathworld.wolfram.com/ClenshawRecurrenceFormula.html
 *
 * Let\verbatim
    S = sum(k = 0..N) c[k] * F[k](x)
    F[n+1](x) = alpha[n](x) * F[n](x) + beta[n](x) * F[n-1](x)
\endverbatim
 * Evaluate <tt>S</tt> with\verbatim
    y[N+2] = y[N+1] = 0
    y[k] = alpha[k] * y[k+1] + beta[k+1] * y[k+2] + c[k]
    S = c[0] * F[0] + y[1] * F[1] + beta[1] * F[0] * y[2]
\endverbatim
 * \e IF <tt>F[0](x) = 1</tt> and <tt>beta(0,x) = 0</tt>, then <tt>F[1](x) =
 * alpha(0,x)</tt> and we can continue the recursion for <tt>y[k]</tt> until
 * <tt>y[0]</tt>, giving\verbatim
    S = y[0]
\endverbatim
 *
 * Evaluating the inner sum\verbatim
 l = n-m; n = l+m
 Sc[m] = sum(l = 0..N-m) C[l+m,m] * q^l * P[l+m,m](t)/P[m,m](t)
 F[l] = q^l * P[l+m,m](t)/P[m,m](t)
\endverbatim
 * Holmes + Featherstone, Eq. (11), give\verbatim
   P[n,m] = sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m))) * t * P[n-1,m] -
            sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3))) * P[n-2,m]
\endverbatim
 * thus\verbatim
   alpha[l] = t * q * sqrt(((2*n+1)*(2*n+3))/
                           ((n-m+1)*(n+m+1)))
   beta[l+1] = - q^2 * sqrt(((n-m+1)*(n+m+1)*(2*n+5))/
                            ((n-m+2)*(n+m+2)*(2*n+1)))
\endverbatim
 * In this case, <tt>F[0] = 1</tt> and <tt>beta[0] = 0</tt>, so the <tt>Sc[m]
 * = y[0]</tt>.
 *
 * Evaluating the outer sum\verbatim
 V = sum(m = 0..N) Sc[m] * q^(m+1) * cos(m*lambda) * P[m,m](t)
   + sum(m = 0..N) Ss[m] * q^(m+1) * cos(m*lambda) * P[m,m](t)
 F[m] = q^(m+1) * cos(m*lambda) * P[m,m](t) [or sin(m*lambda)]
\endverbatim
 * Holmes + Featherstone, Eq. (13), give\verbatim
   P[m,m] = u * sqrt((2*m+1)/((m>1?2:1)*m)) * P[m-1,m-1]
\endverbatim
 * also, we have\verbatim
   cos((m+1)*lambda) = 2*cos(lambda)*cos(m*lambda) - cos((m-1)*lambda)
\endverbatim
 * thus\verbatim
   alpha[m] = 2*cos(lambda) * sqrt((2*m+3)/(2*(m+1))) * u * q
            =   cos(lambda) * sqrt( 2*(2*m+3)/(m+1) ) * u * q
   beta[m+1] = -sqrt((2*m+3)*(2*m+5)/(4*(m+1)*(m+2))) * u^2 * q^2
               * (m == 0 ? sqrt(2) : 1)
\endverbatim
 * Thus\verbatim
 F[0] = q                                [or 0]
 F[1] = cos(lambda) * sqrt(3) * u * q^2  [or sin(lambda)]
 beta[1] = - sqrt(15/4) * u^2 * q^2
\endverbatim
 *
 * Here is how the various components of the gradient are computed
 *
 * Differentiate wrt <tt>r</tt>\verbatim
   d q^(n+1) / dr = (-1/r) * (n+1) * q^(n+1)
\endverbatim
 * so multiply <tt>C[n,m]</tt> by <tt>n+1</tt> in inner sum and multiply the
 * sum by <tt>-1/r</tt>.
 *
 * Differentiate wrt <tt>lambda</tt>\verbatim
   d cos(m*lambda) = -m * sin(m*lambda)
   d sin(m*lambda) =  m * cos(m*lambda)
\endverbatim
 * so multiply terms by <tt>m</tt> in outer sum and swap sine and cosine
 * variables.
 *
 * Differentiate wrt <tt>theta</tt>\verbatim
  dV/dtheta = V' = -u * dV/dt = -u * V'
\endverbatim
 * here <tt>'</tt> denotes differentiation wrt <tt>theta</tt>.\verbatim
   d/dtheta (Sc[m] * P[m,m](t)) = Sc'[m] * P[m,m](t) + Sc[m] * P'[m,m](t)
\endverbatim
 * Now <tt>P[m,m](t) = const * u^m</tt>, so <tt>P'[m,m](t) = m * t/u *
 * P[m,m](t)</tt>, thus\verbatim
   d/dtheta (Sc[m] * P[m,m](t)) = (Sc'[m] + m * t/u * Sc[m]) * P[m,m](t)
\endverbatim
 * Clenshaw recursion for <tt>Sc[m]</tt> reads\verbatim
    y[k] = alpha[k] * y[k+1] + beta[k+1] * y[k+2] + c[k]
\endverbatim
 * Substituting <tt>alpha[k] = const * t</tt>, <tt>alpha'[k] = -u/t *
 * alpha[k]</tt>, <tt>beta'[k] = c'[k] = 0</tt> gives\verbatim
    y'[k] = alpha[k] * y'[k+1] + beta[k+1] * y'[k+2] - u/t * alpha[k] * y[k+1]
\endverbatim
 *
 * Finally, given the derivatives of <tt>V</tt>, we can compute the components
 * of the gradient in spherical coordinates and transform the result into
 * cartesian coordinates.
 **********************************************************************/
#ifndef GEOGRAV_HPP
#define GEOGRAV_HPP

#include <cmath>
#include <limits>
#include <algorithm>
#include <assert.h>
#include <initializer_list>

namespace geograv
{

typedef float real_t;

/**
 * The size of the coefficient vector for the cosine terms.

 Using code copied from GeographicLib https://geographiclib.sourceforge.io/
        and slightly modified to use only static memory.
 *
 * @param[in] N the maximum degree.
 * @param[in] M the maximum order.
 * @return the size of the vector of cosine terms as stored in column
 *   major order.
 **********************************************************************/
constexpr int Csize(int N, int M)
{ return (M + 1) * (2 * N - M + 2) / 2; }

/**
 * The size of the coefficient vector for the sine terms.

 Using code copied from GeographicLib https://geographiclib.sourceforge.io/
        and slightly modified to use only static memory.
 *
 * @param[in] N the maximum degree.
 * @param[in] M the maximum order.
 * @return the size of the vector of cosine terms as stored in column
 *   major order.
 **********************************************************************/
constexpr int Ssize(int N, int M)
{ return Csize(N, M);}// - (N + 1); }

/** Struct to hold the sherical harmonic coefficients
Using code copied from GeographicLib https://geographiclib.sourceforge.io/
       and slightly modified to use only static memory.

       NMAX(int greater than 1): maximum degree and order
       */
template<int NMAX>//NMAX maximum degree and order
struct Coeff{
    double earth_radius;
    double earth_gravity_constant;
    double J2;
    real_t _Cnm[Csize(NMAX, NMAX)];
    real_t _Snm[Ssize(NMAX, NMAX)];
    real_t _sqrttable[std::max(2*NMAX+ 5, 15) + 1];

    /**
    returns the maximum degree and order of the model.
    */
    constexpr int max_degree() const
    { return NMAX; }

    /**
    returns the sqrt of an integer range 0 to std::max(2*NMAX+ 5, 15).
    */
    inline real_t intsqrt(int i) const
    {return _sqrttable[i];}

    /**
     * The one-dimensional index into \e C and \e S.
     *
     * @param[in] n the degree.
     * @param[in] m the order.
     * @return the one-dimensional index.
     **********************************************************************/
    constexpr int index(int n, int m) const
    { return ((((NMAX<<1) - m + 1)*m)>>1) + n; }

    /**
     * An element of \e C.
     *
     * @param[in] k the one-dimensional index.
     * @return the value of the \e C coefficient.
     **********************************************************************/
    real_t Cv(int k) const { return _Cnm[k]; }

    /**
     * An element of \e S.
     *
     * @param[in] k the one-dimensional index.
     * @return the value of the \e S coefficient.
     **********************************************************************/
    real_t Sv(int k) const { return _Snm[k]; }// - (NMAX + 1)]; }


    /**
    cast the Coeff to a Coeff with different maximum degree and order
    */
    template<int OTHER_NMAX>//new Coeff maximum degree and order
    constexpr operator Coeff<OTHER_NMAX>() const{
      //static_assert(OTHER_NMAX<=NMAX,"new max order and degree must be equal or lower");
      Coeff<OTHER_NMAX> other={0,0,0,{0},{0},{0}};
      other.earth_radius= earth_radius;
      other.earth_gravity_constant= earth_gravity_constant;
      other.J2= J2;
      for (int n=0; n<= OTHER_NMAX; n++){
        for (int m=0; m<= n; m++){
          if ((m>NMAX) || (n>NMAX)){
            other._Cnm[other.index(n,m)]=0;
            other._Snm[other.index(n,m)]=0;
          }else {
            other._Cnm[other.index(n,m)]=_Cnm[index(n,m)];
            other._Snm[other.index(n,m)]=_Snm[index(n,m)];
          }
        }
      }
      for (int i=0; i< (std::max(2*OTHER_NMAX+ 5, 15) + 1); i++){
        if (i>std::max(2*NMAX+ 5, 15)){
          other._sqrttable[i]=std::sqrt(double(i));
        }else
          other._sqrttable[i]=_sqrttable[i];
      }
      return other;
    }

};





typedef struct {
    double x;
    double y;
    double z;
} Vector;

// An internal scaling of the coefficients to avoid overflow in
// intermediate calculations.
// copied from GeographicLib https://geographiclib.sourceforge.io/
const real_t scale= real_t(std::pow(real_t(std::numeric_limits<real_t>::radix),
                        -3 * (std::numeric_limits<real_t>::max_exponent < (1<<14) ?
                        std::numeric_limits<real_t>::max_exponent : (1<<14))
                        / 5));

// Move latitudes near the pole off the axis by this amount.
// copied from GeographicLib https://geographiclib.sourceforge.io/
const real_t eps= std::numeric_limits<real_t>::epsilon() *
                      std::sqrt(std::numeric_limits<real_t>::epsilon());




/** Return the gravity potential in International Terrestrial Reference System coordinates, units J/kg.
 @param[in] position_itrs(Above the surface of earth): The location where the gravity is calculated, units m.
 @param[in]  c(): gravity model to use.
 @param[in]  add_pointmass_gravity(): if true, calculate the full gravity, if false, only calculate the non point mass part.
 @param[out] acceleration: The acceleration due to gravity, units m/s^2.

    Using code copied from https://geographiclib.sourceforge.io/ and slightly modified to run in mixed precision and use static memory.
    J2 and J0(point mass) terms are calculated in double precision, while the higher order terms are calculated in single precision.
 */
template<int NMAX>
inline double GeoGrav(Vector position_itrs, Vector& acceleration, const Coeff<NMAX>& c,bool add_pointmass_gravity){
    int N= NMAX;
    int M= NMAX;
    double a= c.earth_radius;
    double x= position_itrs.x;
    double y= position_itrs.y;
    double z= position_itrs.z;
    double
      p= std::hypot(x,y),
      cl = p != 0 ? x / p : 1,  // cos(lambda); at pole, pick lambda = 0
      sl = p != 0 ? y / p : 0,  // sin(lambda)
      r = std::hypot(p,z),
      t = r != 0 ? z / r : 0,   // cos(theta); at origin, pick theta = pi/2
      u = r != 0 ? std::max(p / r, double(eps)) : 1, // sin(theta); but avoid the pole
      q = a / r;
    double
      q2 = q*q,
      uq = u * q,
      uq2 = uq*uq,
      tu = t / u;
    // Initialize outer sum
    real_t vc  = 0, vc2  = 0, vs  = 0, vs2  = 0;   // v [N + 1], v [N + 2]
    // vr, vt, vl and similar w variable accumulate the sums for the
    // derivatives wrt r, theta, and lambda, respectively.
    real_t vrc = 0, vrc2 = 0, vrs = 0, vrs2 = 0;   // vr[N + 1], vr[N + 2]
    real_t vtc = 0, vtc2 = 0, vts = 0, vts2 = 0;   // vt[N + 1], vt[N + 2]
    real_t vlc = 0, vlc2 = 0, vls = 0, vls2 = 0;   // vl[N + 1], vl[N + 2]

    int k;
    for (int m = M; m >= 0; --m) {   // m = M .. 0
      // Initialize inner sum
      real_t
        wc  = 0, wc2  = 0, ws  = 0, ws2  = 0, // w [N - m + 1], w [N - m + 2]
        wrc = 0, wrc2 = 0, wrs = 0, wrs2 = 0, // wr[N - m + 1], wr[N - m + 2]
        wtc = 0, wtc2 = 0, wts = 0, wts2 = 0; // wt[N - m + 1], wt[N - m + 2]
      k = c.index(N, m) + 1;

      for (int n = N; n >= m; --n) {             // n = N .. m; l = N - m .. 0
        real_t w, A, Ax, B, R;    // alpha[l], beta[l + 1]
          w = c.intsqrt((2 * n + 1)) / (c.intsqrt((n - m + 1)) * c.intsqrt((n + m + 1)));
          Ax = real_t(q) * w * c.intsqrt((2 * n + 3));
          A = real_t(t) * Ax;
          B = - real_t(q2) * c.intsqrt((2 * n + 5)) /
            (w * c.intsqrt((n - m + 2)) * c.intsqrt((n + m + 2)));
        R = c.Cv(--k);
        R *= scale;
        w = A * wc + B * wc2 + R; wc2 = wc; wc = w;
          w = A * wrc + B * wrc2 + (n + 1) * R; wrc2 = wrc; wrc = w;
          w = A * wtc + B * wtc2 -  real_t(u)*Ax * wc2; wtc2 = wtc; wtc = w;
        if (m) {
          R = c.Sv(k);
          R *= scale;
          w = A * ws + B * ws2 + R; ws2 = ws; ws = w;
            w = A * wrs + B * wrs2 + (n + 1) * R; wrs2 = wrs; wrs = w;
            w = A * wts + B * wts2 -  real_t(u)*Ax * ws2; wts2 = wts; wts = w;
        }
      }
      // Now Sc[m] = wc, Ss[m] = ws
      // Sc'[m] = wtc, Ss'[m] = wtc
      if (m) {
        real_t v, A, B;           // alpha[m], beta[m + 1]
          v = c.intsqrt((2)) * c.intsqrt((2 * m + 3)) / c.intsqrt((m + 1));
          A = real_t(cl) * v * real_t(uq);
          B = - v * c.intsqrt((2 * m + 5)) / (c.intsqrt((8)) * c.intsqrt((m + 2))) * uq2;
        v = A * vc  + B * vc2  +  wc ; vc2  = vc ; vc  = v;
        v = A * vs  + B * vs2  +  ws ; vs2  = vs ; vs  = v;
          // Include the terms Sc[m] * P'[m,m](t) and Ss[m] * P'[m,m](t)
          wtc += m * real_t(tu) * wc; wts += m * real_t(tu) * ws;
          v = A * vrc + B * vrc2 +  wrc; vrc2 = vrc; vrc = v;
          v = A * vrs + B * vrs2 +  wrs; vrs2 = vrs; vrs = v;
          v = A * vtc + B * vtc2 +  wtc; vtc2 = vtc; vtc = v;
          v = A * vts + B * vts2 +  wts; vts2 = vts; vts = v;
          v = A * vlc + B * vlc2 + m*ws; vlc2 = vlc; vlc = v;
          v = A * vls + B * vls2 - m*wc; vls2 = vls; vls = v;
      } else {
        real_t A, B, qs;
          A = c.intsqrt((3)) * real_t(uq);       // F[1]/(q*cl) or F[1]/(q*sl)
          B = - c.intsqrt((15))/real_t(2) * real_t(uq2); // beta[1]/q
        qs = real_t(q / scale);
        vc = qs * (wc + A * (real_t(cl) * vc + real_t(sl) * vs ) + B * vc2);
          qs /= real_t(r);
          // The components of the gradient in spherical coordinates are
          // r: dV/dr
          // theta: 1/r * dV/dtheta
          // lambda: 1/(r*u) * dV/dlambda
          vrc =   - qs * (wrc + A * (cl * vrc + sl * vrs) + B * vrc2);
          vtc =     qs * (wtc + A * (cl * vtc + sl * vts) + B * vtc2);
          vlc = qs / u * (      A * (cl * vlc + sl * vls) + B * vlc2);
      }
    }

      // Rotate into cartesian (geocentric) coordinates
      double f= c.earth_gravity_constant/c.earth_radius;
      acceleration.x = cl * (u * vrc + t * vtc) - sl * vlc;
      acceleration.y = sl * (u * vrc + t * vtc) + cl * vlc;
      acceleration.z =       t * vrc - u * vtc            ;
      //Calculate J2 term in double precission
      double t2= t*t;
      double q4= q2*q2;
      double temp= (5.0L*(t2)-1.0L)*(3.0L/2.0L)*q4;
      double v31= temp*(x/r);
      double w31= temp*(y/r);
      double v20= q2*q*0.5L*(3*t2-1);
      double v30= ((5.0L/3.0L)*(t2)-1.0L)*t*(3.0L/2.0L)*q4;
      double J2_unnorm= std::sqrt(5.0L)*c.J2;
      double j2x= -J2_unnorm*v31*(1.0L/a);
      double j2y= -J2_unnorm*w31*(1.0L/a);
      double j2z= -J2_unnorm*v30*(3.0L/a);
      double j2pot= J2_unnorm*v20;
      acceleration.x +=j2x;
      acceleration.y +=j2y;
      acceleration.z +=j2z;
      double pot= j2pot+double(vc);

      if(add_pointmass_gravity){
        acceleration.x += -(x/(r*r))*q;
        acceleration.y += -(y/(r*r))*q;
        acceleration.z += -(z/(r*r))*q;
        pot += q;
      }

      acceleration.x*=f;
      acceleration.y*=f;
      acceleration.z*=f;
    return pot*f;
  }

}
#endif /* GEOGRAV_HPP */
