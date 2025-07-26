
/*
* * 
* * Filename: Vibrato LFO
* *
* * Description:
* *   A low-frequency oscillator (LFO) for creating vibrato effects.
* *   This LFO can be used to modulate parameters such as pitch or amplitude
* *   to create a vibrato effect in audio synthesis. It works in combination
* *   with a Farrow intepolator to achieve smooth modulation.
* *
* *
* * Author:
* * JEP J. Enrique Peraza
*/
#pragma once
#include <vector>
#include <tuple>
#include <cmath>
#include <iostream>

namespace sig::wg
{
  /* ------------------------------------------------------------ */
/*  Shared helper: vibrato LFO                                   */
/* ------------------------------------------------------------ */
struct VibLFO {
    double depth{0.0};                  // The depth of the vibrato
    double rate{0.0};                   // The vibrato rate in Hz.
    double phase{0.0};                  // The current phase of the LFO
    // Enable vibratio modulation.
    inline void Enable(                 // Enable LFO modulation.
      double d,                         // The depth of the vibrato
      double r) noexcept                // The rate of the vibrato in Hz.
    {                                   // ~~~~~~~~~ Enable ~~~~~~~~~~ //
      depth=d;                          // Set the depth of the vibrato
      rate=r;                           // Set the rate of the vibrato in Hz.
      phase=0.0;                        // Reset the phase of the LFO
    }                                   // ~~~~~~~~~ Enable ~~~~~~~~~~ //
    inline void Disable(void){depth=0.0;}// Disable vibrato modulation.
    inline double Tick(double fs)       // Tick the LFO we have to update
    {                                   // ~~~~~~~~~ Tick ~~~~~~~~~~ //
      if(depth<=0||rate<=0) return 0.0; // If no modulation parameters return 0.0
      // Update the phase of the LFO    // Calculate the new phase
      phase=std::fmod(phase+6.28318530718*rate/fs,6.28318530718);
      // Return the current value of the LFO
      return depth*std::sin(phase);     // Oscillate using a sine wave
    }                                   // ~~~~~~~~~ Tick ~~~~~~~~~~ //
};
}
