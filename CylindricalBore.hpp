 /* 
 * *
 * * Filename: CylindricalBore.hpp
 * * Description:
 * * This file contains a CylindricalBore class that represents a cylindrical bore
 * * in a waveguide network. It can make clarinets and on, not very good with voices.
 * * It is used to model the acoustic properties of a cylindrical bore as a waveguide.
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <cassert>
#include <algorithm>
#include <cmath>
#include "VibLFO.hpp"
#include "DelayBranch.hpp"
#include "WGFilters.hpp"
#include "FilterFactory.hpp"
namespace sig::wg
{
  template <typename T=double,
  size_t MaxLen=1<<15,
  size_t K=3,
  size_t P=3>
  class CylindricalBore: public Node
  {
    using Branch=DelayBranch<T,MaxLen,K,P>;
    public:
    bool Prepare(                       // Prime the conical bore
      double fs,                      // Sample rate in Hz.
      double st=0,                  // Stages of dispersion.
      double D=48.0,                  // Integer delay in samples.
      double mt=0.0,                  // Thiran fractional delay in samples.
      double mf=0.0,                  // Farrow mod fractional delay
      size_t o=K) noexcept            // Interpolator bank order
    {                                   // ~~~~~~~~~ Prepare ~~~~~~~~~~~~~~~~~ //
      if (fs<=0.0||fs>=0.5*fs||st<0.0) return false;// Sanitize input
      if (D>0.0&&D<K) return false;     // D must be larger than interpolation order.
      this->fs=fs;                      // Set the sample rate
      this->st=st;                      // Set the reed stiffness coefficient.
      this->D=D;                        // Set the integer delay in samples.
      mut=mt;                           // Set the fractional delay for Thiran.
      muf=mf;                           // Set the fractional delay for Farrow.
      order=o;                          // Set the order of the interpolation filters.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Prepare the loss filter chain
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      chain.Prepare(fs,fs*0.25);        // Prepare the loss filter chain
      size_t L=D?static_cast<size_t>(D):static_cast<size_t>(fs/(2.0*f0)+0.5);
      up.Prepare(L,mut,muf);           // Prepare the upstream delay branch
      dn.Prepare(L,mut,muf);           // Prepare the downstream delay branch
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Prime WG according to Group Delay of these lengths
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      size_t gd=up.GroupDelay(L,o,o);// Get the group delay of the mouthpiece
      for (size_t i=0;i<gd;++i)         // For the length of the Group Delay...
      {                                 // Prime the waveguide branches....
        up.Write(T(0));                 // Write zero to the upstream delay branch
        up.Propagate(1);                // Propagate the zero sample through the upstream delay branch
        dn.Write(T(0));                 // Write zero to the downstream delay branch
        dn.Propagate(1);                // Propagate the zero sample through the downstream delay branch
      }                                 // Done priming waveguide branches
      prepared=true;                    // A loaded gun.
      return true;                      // Preparation successful
    }                                   // ~~~~~~~~~ Prepare ~~~~~~~~~~~~~~~~~ //
    inline void Blow(T p) noexcept      // Deposit a sample into whe WG
    {                                   // ~~~~~~~~~ Blow ~~~~~~~~~~~~~~~~~~~~~~ //
        mth=p;                          // Load the mouth pressure.
    }                                   // ~~~~~~~~~ Blow ~~~~~~~~~~~~~~~~~~~~~~ //
    inline void Propagate(size_t n) noexcept      // Process the waveguide
    {                                   // ~~~~~~~~~ Proapagate ~~~~~~~~~~~~~~~~~ //
      for (size_t i=0;i<n;++i)          // For # of samples to propagate...
      {                                 // Propagate the waveguide branches
        double doff=lfo.Tick(fs);       // LFO phase offset
        up.SetFractionalDelay(mut+doff); // Set the fractional delay for Thiran with vibrato modulation
        dn.SetFractionalDelay(mut+doff); // Set the fractional delay for Thiran with vibrato modulation
        up.SetMuFarrow(muf+doff);       // Set the fractional delay for Farrow with vibrato modulation
        dn.SetMuFarrow(muf+doff);
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Model: Simple linear reed and radiation junction.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //         
        double pup=up.Read();           // Read the upstream wave samples
        double pdn=dn.Read();           // Read the downstream wave samples
        double dp=mth-pup;              // Pressure wave differential
        double rfl=(dp>0.0?dp*dp*invk:0);// Volumetric flow wave from reed.
        double preed=mth-rfl;          // Pressure wave at reed junction
        double prad=-orfl*pdn;          // Radiated pressure wave.
        // Output waves to write to branches.
        double outu=chain.Process(preed+prad);// Process the upstream wave with loss
        double outd=chain.Process(-prad);// Process the radiated wave.
        // Propagate the waves through the branches.
        up.Write(outu);                // Write the upstream wave sample
        up.Propagate(1);               // Propagate the upstream wave sample
        dn.Write(outd);                // Write the downstream wave sample
        dn.Propagate(1);               // Propagate the downstream wave sample
      }                                // Done propagating the waveguide branches
    }                                   // ~~~~~~~~~ Propagate ~~~~~~~~~~~~~~~~~ //
    // Output API:
    inline T Pressure(void) const noexcept // Get the current pressure output
    {                                   // ~~~~~~~~~ Pressure ~~~~~~~~~~~~~~~~~ //
      return static_cast<T>(mth);       // Return the mouth pressure as output
    }                                   // ~~~~~~~~~ Pressure ~~~~~~~~~~~~~~~~~ //
    // Getters and Setters API:
    inline double GetSampleRate(void) const noexcept
    {                                   // ~~~~~~~~~ GetSampleRate ~~~~~~~~~~~~~~~~~ //
      return fs;                        // Return the sample rate in Hz.
    }                                   // ~~~~~~~~~ GetSampleRate ~~~~~~~~~~~~~~~~~ //
    inline void SetSampleRate(double s) noexcept
    {                                   // ~~~~~~~~~ SetSampleRate ~~~~~~~~~~~~~~~~~ //
      fs=s;                              // Set the sample rate in Hz.
      chain.Prepare(s,fs*0.25);          // Prepare the loss filter chain with the new sample rate
    }                                   // ~~~~~~~~~ SetSampleRate ~~~~~~~~~~~~~~~~~ //
    inline void SetSamplingFrequency(double fs) noexcept // Alias for SetSampleRate
    {                                   // ~~~~~~~~~ SetSamplingFrequency ~~~~~~~~~ //
      SetSampleRate(fs);                // Call SetSampleRate
    }                                   // ~~~~~~~~~ SetSamplingFrequency ~~~~~~~~~ //
    inline void SetCutoffFrequency(double fc) noexcept // Set loss filter cutoff
    {                                   // ~~~~~~~~~ SetCutoffFrequency ~~~~~~~~~~~ //
      chain.Prepare(fs, fc);            // Re-prepare the loss filter with new cutoff
    }                                   // ~~~~~~~~~ SetCutoffFrequency ~~~~~~~~~~~ //
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Vibrato modulation API:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    inline void EnableVibrato(
      double d,                        // The depth of the vibrato
      double r) noexcept               // The rate of the vibrato in Hz.
    {                                  // ~~~~~~~~~ EnableVibrato ~~~~~~~~~~~~~~~ //
      lfo.Enable(d,r);                 // Enable the vibrato LFO with the given depth and rate.
    }                                  // ~~~~~~~~~ EnableVibrato ~~~~~~~~~~~~~~~ //
    inline void DisableVibrato(void) {lfo.Disable();}// Disable vib modulation.

    // Set Delay of the All-Pass approximator
    inline void SetFractionalDelay(double mt) noexcept
    {                                   // ~~~~~~~~~ SetFractionalDelay ~~~~~~~~~~ //
      mut=mt;                           // Set the fractional delay for Thiran.
      up.SetFractionalDelay(mut);       // Set the fractional delay for Thiran.
      dn.SetFractionalDelay(mut);       // Set the fractional delay for Thiran.
    }                                   // ~~~~~~~~~ SetFractionalDelay ~~~~~~~~~~ //
    inline void SetMuFarrow(double mf) noexcept
    {                                   // ~~~~~~~~~ SetMuFarrow ~~~~~~~~~~~~~~~ //
      muf=mf;                           // Set the fractional delay for Farrow.
      up.SetMuFarrow(muf);              // Set the fractional delay for Farrow.
      dn.SetMuFarrow(muf);              // Set the fractional delay for Farrow.
    }                                   // ~~~~~~~~~ SetMuFarrow ~~~~~~~~~~~~~~~ //
    inline double GetMuFarrow(void) const noexcept
    {                                   // ~~~~~~~~~ GetMuFarrow ~~~~~~~~~~~~~~~ //
      return muf;                       // Return the fractional delay for Farrow.
    }                                   // ~~~~~~~~~ GetMuFarrow ~~~~~~~~~~~~~~~ //
    inline void SetOrder(size_t o) noexcept
    {
      order=o;                          // Set the order of the Farrow deinterpolator.
    }                                   // ~~~~~~~~~ SetOrder ~~~~~~~~~~~~~~~~~ //
    inline size_t GetOrder(void) const noexcept
    {                                   // ~~~~~~~~~ GetOrder ~~~~~~~~~~~~~~~~~ //
      return order;                     // Return the order of the Farrow deinterpolator.
    }                                   // ~~~~~~~~~ GetOrder ~~~~~~~~~~~~~~~~~ // 
    inline void SetStages(double s) noexcept
    {                                   // ~~~~~~~~~ SetStages ~~~~~~~~~~~~~~~~~ //
      st=s;                              // Set the stages of dispersion.
    }                                   // ~~~~~~~~~ SetStages ~~~~~~~~~~~~~~~~~ //
    inline double GetStages(void) const noexcept
    {                                   // ~~~~~~~~~ GetStages ~~~~~~~~~~~~~~~~~ //
      return st;                        // Return the stages of dispersion.
    }                                   // ~~~~~~~~~ GetStages ~~~~~~~~~~~~~~~~~ //
    inline void PitchBend(double c) noexcept 
    {                                   // ~~~~~~~~~ PitchBend ~~~~~~~~~~~~~~~~~ //
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ // Using the pitch bend in cents, 
      // Get the -3db cutoff frequency for the pitch bend.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
       double fc=std::pow(2.0,c/1200.0);// Calculate the pitch bend factor
       Prepare(fs,st,f0*fc,mut,muf,order);// Re-prepare the conical bore with the new pitch bend factor
    }                                    // Pitch-bend 
    inline void Clear(void) noexcept     // Clear the waveguide state
    {                                   // ~~~~~~~~~ Clear ~~~~~~~~~~~~~~~~~~~~~~~~ //
      up.Clear();                       // Clear the upstream delay branch
      dn.Clear();                       // Clear the downstream delay branch
      chain.Clear();                    // Clear the loss filter chain
      lfo.Disable();                    // Disable the vibrato LFO
      prepared=false;                   // Reset the prepared flag
      mth=0.0;                          // Reset the mouth pressure
      fs=48000.0;                       // Reset the sample rate
      f0=440.0;                         // Reset the fundamental frequency
      D=0.0;                            // Reset the integer delay in samples
      mut=0.0;                          // Reset the fractional delay in samples
      muf=0.0;                          // Reset the fractional delay for Farrow
      order=K;                          // Reset the order of the Farrow deinterpolator
      st=5.0;                           // Reset the stages of dispersion
      pos=0.15;                         // Reset the position
    }                                   // ~~~~~~~~~ Clear ~~~~~~~~~~~~~~~~~~~~~~~~ //
    inline void SetPosition(double p) noexcept
    {
      pos=p;                            // Set the position in the string
    }
    inline double GetPosition(void) const noexcept
    {
      return pos;                       // Return the position in the string
    }
    // Members
    private:
      Branch up,dn;                     // Upstream and downstream delay branches
      LossChain<T> chain;               // Loss filter chain
      VibLFO lfo;                       // Vibrato LFO for modulation
      double fs{48000.0},               // Sample rate in Hz
        f0{440.0},                      // Fundamental frequency in Hz
        mth{0.0},                       // Mouth pressure
        orfl{0.9},                      // Open wave impedance mismatch reflection coefficient.
        // Admittance coefficient as a function of the reflection 
        // of the wave due to an opening.
        invk{1.0/9.0},                  // Inverse of the open wave impedance mismatch reflection coefficient.
        st{5.0},                        // Reed stiffness coefficient
        D{0.0},                         // Integer delay in samples
        alpha{0.18},                    // Dispersion coefficient (0.0 to 0.5)
        mut{0.0},                       // Thiran fractional delay in samples
        muf{0.0},                       // Farrow mod fractional delay
        pos{0.15},                      // Position in the string
        order{K};                       // Interpolator bank order
      bool prepared{false};             // Preparation status
  };
}
