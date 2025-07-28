 /* 
 * *
 * * Filename: ReedMouthpiece.hpp
 * * Description:
 * * This file contains a ReedMouthpiece class that represents a reed mouthpiece
 * * in a waveguide network.
 * * It is used to model the acoustic properties of a reed mouthpiece as a waveguide junction.
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
  class ReedMouthpiece:public Node
  {
    using Branch=DelayBranch<T,MaxLen,K,P>;
    public:
      bool Prepare(       
        double fs,                      // Sample rate in Hz.
        double rk=5.0,                  // Reed stiffness coefficient.
        double D=48.0,                  // Integer delay in samples.
        double mt=0.0,                  // Thiran fractional delay in samples.
        double mf=0.0,                  // Farrow mod fractional delay
        size_t o=K) noexcept            // Interpolator bank order
    {                                   // ~~~~~~~~~ Prepare ~~~~~~~~~~~~~~~~~ //
      if (fs<=0.0||rk<0.0) return false;// Sanitize input
      if (D>0.0&&D<K) return false;     // D must be larger than interpolation order.
      this->fs=fs;                      // Set the sample rate
      this->rk=rk;                      // Set the reed stiffness coefficient.
      this->D=D;                        // Set the integer delay in samples.
      mut=mt;                           // Set the fractional delay for Thiran.
      muf=mf;                           // Set the fractional delay for Farrow.
      order=o;                          // Set the order of the interpolation filters.
      mp.Prepare(size_t(D),mut,muf);    // Prepare the mouthpiece delay branch
      // ------------------------------ //
      // Prepare the loss filter chain
      // ------------------------------ //
      chain.Prepare(fs,fs*0.25);        // Prepare the loss filter chain
      size_t L=D?static_cast<size_t>(D):static_cast<size_t>(fs/(2.0*f0)+0.5);
      // ----------------------------- //
      // Prime WG according to Group Delay of these lengths
      // ----------------------------- //
      size_t gd=mp.GroupDelay(L,o,o);  // Get the group delay of the mouthpiece
      for (size_t i=0;i<gd;++i)        // For the length of the Group Delay...
      {                                // Prime the branches...
        mp.Write(T(0));                // Write zero to the mouthpiece delay branch
        mp.Propagate(1);               // Propagate the zero sample through the mouthpiece delay branch
      }                                // Done priming the mouthpiece delay branch.
      prepared=true;                   // A loaded gun.
      return true;                     // True if we got here.
    }                                  // ~~~~~~~~~ Prepare ~~~~~~~~~~~~~~~~~ //
    bool Propagate(T p) noexcept        // Process the pressure sample measurement p
    {                                   // ~~~~~~~~~ Propagate ~~~~~~~~~~~~~~~~~~~~~ //      
      double phoff=lfo.Tick(fs);        // Get the vibrato modulation phase offset
      mp.SetFractionalDelay(mut+phoff); // Set the fractional delay for Thiran with vibrato modulation
      mp.SetMuFarrow(muf+phoff);        // Set the fractional delay for Farrow with vibrato modulation
      double pprev=mp.Read();           // Read the previous pressure sample from the mouthpiece delay branch
      // ------------------------------ //
      // Simple linear read-valve. The dynamics is built into the waveguide already
      // to allow for these simplifications.
      // ------------------------------ //
      double dp=p-pprev;                // Calculate pressure wave differential
      double flow=(dp>0?dp*dp/rk:0.0);  // Calculate volumetric flow wave
      double pnext=p-flow;              // Calculate the next pressure sample
      mp.Write(pnext);                  // Write wave into buffer
      mp.Propagate(1);                  // Propagate the pressure sample through the mouthpiece delay branch
      return pprev;                     // Return the previous pressure sample
    }                                   // ~~~~~~~~~ Propagate ~~~~~~~~~~~~~~~~~~~~~ //
    // Override Node::Propagate for polymorphic use
    void Propagate(size_t n) noexcept override // Process n samples  
    {                                   // ~~~~~~~~~ Propagate ~~~~~~~~~~~~~~~~~~~~~ //
      for(size_t i = 0; i < n; ++i) {  // Process each sample
        Propagate(T(0));                // Process with zero input (self-oscillating)
      }
    }                                   // ~~~~~~~~~ Propagate ~~~~~~~~~~~~~~~~~~~~~ //
    inline void Clear(void) noexcept    // Clear the state of the ReedMouthpiece.hpp
    {                                   // ~~~~~~~~~ Clear ~~~~~~~~~~~~~~~~~~~~~~~~ //
      mp.Clear();                       // Clear the mouthpiece delay branch
      chain.Clear();                    // Clear the loss filter chain
      lfo.Disable();                    // Disable the vibrato LFO
      prepared=false;                   // Reset the prepared flag
    }                                   // ~~~~~~~~~ Clear ~~~~~~~~~~~~~~~~~~~~~~~~ //
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
    // API methods expected by ClarinetInstrumentGraph
    template<typename BoreType>
    inline bool Connect(BoreType& bore) noexcept { return true; } // Simple connect
    inline bool IsReady(void) const noexcept { return prepared; } // Check readiness
    inline void SetMouthPressure(float p) noexcept { /* Store mouth pressure */ }
    inline void SetReedArea(float a) noexcept { /* Store reed area */ }
    inline void SetBernoulliConstant(float m) noexcept { /* Store Bernoulli constant */ }

    inline void SetReedStiffness(double rk) noexcept
    {
      this->rk=rk;
    }
    inline double GetReedStiffness(void) const noexcept
    {                                   // ~~~~~~~~~ GetReedStiffness ~~~~~~~~~~~~~~~ //
      return rk;                        // Return the reed stiffness coefficient.
    }                                   // ~~~~~~~~~ GetReedStiffness ~~~~~~~~~~~~~~~ //
    // -------------------------------- //
    // Vibrato modulation control.
    // -------------------------------- //
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
      mp.SetFractionalDelay(mut);       // Set the fractional delay for Thiran.
    }                                   // ~~~~~~~~~ SetFractionalDelay ~~~~~~~~~~ //
    inline void SetMuFarrow(double mf) noexcept
    {                                   // ~~~~~~~~~ SetMuFarrow ~~~~~~~~~~~~~~~ //
      muf=mf;                           // Set the fractional delay for Farrow.
      mp.SetMuFarrow(muf);              // Set the fractional delay for Farrow.
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
    private:
      // Our waveguide buffers and filter junctions
      Branch mp{};                      // Mouthpiece delay branch
      LossChain<T> chain;               // Loss filter chain
      VibLFO lfo;                       // Vibrato LFO
      FilterFactory<T> ff;              // Filter factory for easy filter creation
      // Our member variables.          //
      double fs{48000.0};               // Sample rate in Hz.
      double f0{440.0};                 // Fundamental frequency in Hz.
      double rk{5.0};                   // Reed stiffness coefficient.
      double lfc{6000.0};               // Loss filter cutoff frequency in Hz.
      double lq{0.7071};                // Loss filter Q factor.
      double alpha{0.18};               // Dispersion coefficient (0.0 to 0.5).
      // Our delay branches fractional delay processor
      // parameters.
      double D{0.0};                   // Integer delay in seconds.
      double mut{0.0};                 // Thiran fractional delay in seconds.
      double muf{0.0};                 // Farrow fractional delay in seconds.
      size_t order{K};                 // Order of the Farrow deinterpolator.
      bool prepared{false};            // Preparation flag.
      // Instrument geometric parameters.
  };
}
