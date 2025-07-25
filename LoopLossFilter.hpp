/*
* *
* * Filename: LoopLossFilter.hpp
* *
* * Description: One-Pole exponential smoother (Low-Pass) used to model frequency dependent energy loss
* *  in waveguide conical, tube and string like structures.
* *
* * y[n] = y[n-1]+alpha(x[n]-y[n])
* *   where alpha is the smoothing factor and is defined as:
* *   alpha = dt (RC + dt) with the time constant RC = 1/2*pi*fc
* *   and dt is the sample period.
* *
* * Author:
* * JEP J. Enrique Peraza
* *
*/
#pragma once
#include <cmath>
namespace sig::wg
{
  template <typename T=float>
  class LoopLossFilter
  {
    public:
      void Prepare(                     // Prepare the loop loss filter
        double fs,                      // Sample rate
        double fc) noexcept             // -3dB cutoff frequency
    {                                   // ~~~~~~~~~ Prepare ~~~~~~~~~ //
        if (fc<=0.0 || fs<=0.0 || fc>=0.5*fs) return;
        const T RC=static_cast<T>(1.0/2.0*M_PI*fc); // Time constant
        const T dt=static_cast<T>(1.0/fs);// Sample period
        alpha=dt/(RC+dt);               // Smoothing factor (0 < alpha < 1)'
        y1=static_cast<T>(0.0);         // Reset the previous output sample
    }                                   // ~~~~~~~~~ Prepare ~~~~~~~~~ //
    inline T ProcessSample(T x) noexcept
    {                                   // --------- ProcessSample --------- //
    y1+=alpha*(x-y1);                   // Apply the smoothing filter
    return y1;                          // Return the output sample
    }                                   // --------- ProcessSample --------- //
    inline void Clear(void) noexcept
    {                                   // ~~~~~~~~~ Clear ~~~~~~~~~ //
      y1=static_cast<T>(0.0);           // Reset the previous output sample
    }                                   // ~~~~~~~~~ Clear ~~~~~~~~~ //
    inline T GetAlpha(void) const noexcept
    {                                   // ~~~~~~~~~ GetAlpha ~~~~~~~~~ //
        return alpha;                   // Return the smoothing factor
    }                                   // ~~~~~~~~~ GetAlpha ~~~~~~~~~ //
    inline T GetOutput(void) const noexcept
    {
        return y1;                       // Return the last output sample
    }                                    // ~~~~~~~~~ GetOutput ~~~~~~~~~ //
    private:
      double alpha{0.0};  // Smoothing factor
      T y1{0.0}; // Previous output sample
  };
} // namespace sig::wg
