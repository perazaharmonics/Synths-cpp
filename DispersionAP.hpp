/*  DispersionAP.hpp
 *  -------------------------------------------
 *  First-order all-pass:  H(z) = (z?¹ - a) / (1 - a z?¹)
 *  |a| < 1  gives a phase slope that approximates
 *  stiffness dispersion in strings or jet dispersion
 *  in wind instruments.  Chain several stages for
 *  stronger inharmonicity.
 *  -------------------------------------------
 */
#pragma once
namespace sig::wg
{
  template<typename T=float>
  class DispersionAllPass
  {
    public:
      // Alpha in the range 0 ... 0.5 for subtle dispersion.
      // negative values flip phase slope (rarely needed).
      void Prepare(T alpha) noexcept
      {
        a=alpha;
        z1=T(0); // Reset the delay line
      }
      inline T ProcessSample(T x) noexcept
      {                            // --------- ProcessSample --------- //
        const T y=-a*x+z1+a*z1;    // Direct Form I
        z1=x;                      // Update the delay line
        return y;                  // Return the output sample
      }                            // --------- ProcessSample --------- //
      void Clear(void) noexcept
      {                            // Clear the state of the all-pass filter
        z1=T(0);                   // Reset the delay line
      }                            // ---------- Clear --------- //
      // Access the dispersion coefficient
      inline void SetCoefficient(T alpha) noexcept { a = alpha; }
      inline T GetCoefficient(void) const noexcept { return a; }
    private:
      T a=T(0);                    // Coefficient of the all-pass filter
      T z1=T(0);                   // z?¹
  };
}
