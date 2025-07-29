
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
    double fs{48000.0};                 // Sample rate in Hz.
    // Enable vibrato modulation.
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
      if (fs<=0||fs>0.5*this->fs) return 0.0; // If sample rate is invalid, return 0.0
      this->fs=fs;                     // Set the sample rate
      // Update the phase of the LFO    // Calculate the new phase
      phase=std::fmod(phase+6.28318530718*rate/fs,6.28318530718);
      // Return the current value of the LFO
      return depth*std::sin(phase);     // Oscillate using a sine wave
    }                                   // ~~~~~~~~~ Tick ~~~~~~~~~~ //
    inline void SetPhase(double p) noexcept
    {                                   // ~~~~~~~~~ SetPhase ~~~~~~~~~~ //
      phase=std::fmod(p,6.28318530718); // Set the phase of the LFO
    }                                   // ~~~~~~~~~ SetPhase ~~~~~~~~~~ //
    inline double GetPhase(void) const noexcept
    {                                   // ~~~~~~~~~ GetPhase ~~~~~~~~~~ //
      return phase;                     // Return the current phase of the LFO
    }                                   // ~~~~~~~~~ GetPhase ~~~~~~~~~~ //
    inline void SetDepth(double d) noexcept
    {                                   // ~~~~~~~~~ SetDepth ~~~~~~~~~~ //
      depth=d;                          // Set the depth of the vibrato
    }                                   // ~~~~~~~~~ SetDepth ~~~~~~~~~~ //
    inline double GetDepth(void) const noexcept
    {                                   // ~~~~~~~~~ GetDepth ~~~~~~~~~~ //
      return depth;                     // Return the current depth of the vibrato
    }                                   // ~~~~~~~~~ GetDepth ~~~~~~~~~~ //
    inline void SetRate(double r) noexcept
    {                                   // ~~~~~~~~~ SetRate ~~~~~~~~~~ //
      rate=r;                           // Set the rate of the vibrato in Hz.
    }                                   // ~~~~~~~~~ SetRate ~~~~~~~~~~ //
    inline double GetRate(void) const noexcept
    {                                   // ~~~~~~~~~ GetRate ~~~~~~~~~~ //
      return rate;                      // Return the current rate of the vibrato
    }                                   // ~~~~~~~~~ GetRate ~~~~~~~~~~ //
    inline double GetSampleRate(void) const noexcept
    {                                   // ~~~~~~~~~ GetSampleRate ~~~~~~~~~~ //
      return fs;                        // Return the current sample rate of the LFO
    }                                   // ~~~~~~~~~ GetSampleRate ~~~~~~~~~~ //
    inline void SetSampleRate(double s) noexcept
    {                                   // ~~~~~~~~~ SetSampleRate ~~~~~~~~~~ //
      if (s<=0.0||s>=0.5*fs) return;    // If sample rate is invalid, do nothing
      fs=s;                             // Set the sample rate
    }                                   // ~~~~~~~~~ SetSampleRate ~~~~~~~~~~ //
    inline void Reset(void) noexcept
    {                                   // ~~~~~~~~~ Reset ~~~~~~~~~~ //
      phase=0.0;                        // Reset the phase of the LFO
      depth=0.0;                        // Reset the depth of the vibrato
      rate=0.0;                         // Reset the rate of the vibrato
    }                                   // ~~~~~~~~~ Reset ~~~~~~~~~~ //
};
} // namespace sig::wg
