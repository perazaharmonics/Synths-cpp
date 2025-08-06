/*
* *
* * Filename: WaveguideOscillator.hpp
* *
* * Description:
* *  Normalised second‑order digital waveguide oscillator.  This
* *  oscillator is implemented as a simple two‑pole recursive
* *  system whose poles lie on the unit circle at the desired
* *  frequency.  The recursion coefficient is derived from the
* *  cosine of the angular frequency and normalised so that the
* *  amplitude of the output sinusoid remains constant across
* *  frequency.  The oscillator has no inputs; call Process() to
* *  generate successive samples.  Frequency and sample rate can
* *  be changed at run‑time via the provided setters.  Initial
* *  conditions are chosen so that the first output sample is
* *  sin(w) when prepared.
* *
* * Author:
* *  JEP  J. Enrique Peraza
* *
*/
#pragma once
#include <cmath>

namespace sig::wg
{
  template <typename T=float>
  class WaveguideOscillator
  {
    public:
      WaveguideOscillator(void) noexcept = default; // Default constructor
      ~WaveguideOscillator(void) noexcept = default; // Default destructor
      bool Prepare(
        double sr,                    // Sample rate (Hz)
        double f0) noexcept           // Oscillation frequency (Hz)
      {                               // ~~~~~~~~~ Prepare ~~~~~~~~~ //
        if (sr<=0.0||f0>0.5*fs) return false; // Sanity check
        fs=sr;                        // Save sample rate
        f=f0;                         // Save frequency
        // Compute recursion coefficient from frequency
        double w=2.0*std::acos(-1.0)*f0/fs; // Angular frequency
        k=T(2)*static_cast<T>(std::cos(w)); // k = 2*cos(w)
        // Initialise state such that the first output equals sin(w)
        y1=T(0);                      // y[n−1] = 0
        y2=-static_cast<T>(std::sin(w)); // y[n−2] = −sin(w)
        return true;                  // Prepared successfully
      }                               // ~~~~~~~~~ Prepare ~~~~~~~~~ //
      inline void SetFrequency(double f0) noexcept
      {                               // ~~~~~~~~~ SetFrequency ~~~~ //
        f=f0;                         // Save new frequency
        double w=2.0*std::acos(-1.0)*f0/fs; // Compute angular frequency
        k=T(2)*static_cast<T>(std::cos(w)); // Update recursion coefficient
      }                               // ~~~~~~~~~ SetFrequency ~~~~ //
      inline void SetSampleRate (double sr) noexcept
      {                                 // ~~~~~~~~~ SetSampleRate ~~~~ //
        fs=sr;                 // Save new sample rate
        // Recompute coefficient at current frequency
        SetFrequency(f);
      }                                 // ~~~~~~~~~ SetSampleRate ~~~~ //
      inline T Process (void) noexcept
      {                                 // ~~~~~~~~~ Process ~~~~~~~~~ //
        // Compute next oscillator sample using the two‑pole recursion
        T y0=k*y1-y2;                   // y[n] = 2*cos(w)*y[n−1] − y[n−2]
        y2=y1;                          // Shift history
        y1=y0;                          // Save current sample
        return y0;                      // Return oscillator output
      }                                 // ~~~~~~~~~ Process ~~~~~~~~~ //
      inline void Clear(void) noexcept
      {                                 // ~~~~~~~~~ Clear ~~~~~~~~~ //
        y1=T(0);                        // Reset state
        y2=T(0);                        // Reset state
      }                                 // ~~~~~~~~~ Clear ~~~~~~~~~ //
      // Accessors
      inline double GetFrequency(void) const noexcept { return freq; }
      inline double GetSampleRate(void) const noexcept { return fs; }
      inline T GetCoefficient(void) const noexcept { return k; }
    private:
      double fs{48000.0};               // Sample rate in Hz
      double f{440.0};               // Oscillation frequency in Hz
      T k{T(0)};                        // Recursion coefficient = 2*cos(w)
      T y1{T(0)};                       // y[n−1]
      T y2{T(0)};                       // y[n−2]          
  };
} // namespace sig::wg