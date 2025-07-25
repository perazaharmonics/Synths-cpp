/*
* *
* * Filename: WGFilters.hpp
* *
* * Description:
* *  This file contains a handy way to access common waveguide filters for the waveguide layer.
* *  It includes:
* *
* *  - Dispersion All-Pass
* *  - Loop Loss Filter
* *  - Radiation Chain
*/
#pragma once
#include <array>
#include <random>
#include "BiQuad.hpp"
#include "FilterFactory.hpp"
#include "DelayBranch.hpp"
#include "DispersionAP.hpp"
#include "LoopLossFilter.hpp"

namespace sig::wg
{
  template <typename T=float>
  class LossChain
  {
    LoopLossFilter<T> lp1,lp2; // Two loop loss filters
    BiQuad<T> shelf;
    FilterFactory<T> ff; // Filter factory for easy filter creation
    public:
      void Prepare(double fs, double fc) noexcept
      {                                 // Prepare the loss chain with sample rate and cutoff frequency
        lp1.Prepare(fs,fc);             // Prepare the first loop loss filter
        lp2.Prepare(fs,fc);             // Prepare the second loop loss filter
        shelf=ff.HighShelf(fs,fs*0.25,-1.0,0.7);// Prepare the high shelf filter
      }                                 // ~~~~~~~~~ Prepare ~~~~~~~~~ //
      inline T Process(T x) noexcept
      {                                 // ~~~~~~~~~ Process ~~~~~~~~~ //
        return shelf.ProcessSample(lp2.ProcessSample(lp1.ProcessSample(x)));
      }                                 // ~~~~~~~~~ Process ~~~~~~~~~ //
      inline void Clear(void) noexcept
      {                                    // ~~~~~~~~~ Clear ~~~~~~~~~ //
          lp1.Clear();                    // Clear the first loop loss filter
          lp2.Clear();                    // Clear the second loop loss filter
      }                                    // ~~~~~~~~~ Clear ~~~~~~~~~ //
  };

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  // piston High-Pass + bell LP simple radiation model
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  template <typename T=float>
  class RadiationChain
  {
    BiQuad<T> hp,lp;                    // High-pass and low-pass filters         
    FilterFactory<T> ff;                // Filter factory for easy filter creation
    public:
      void Prepare(double fs, double fc) noexcept
      {                                 // ~~~~~~~~~ Prepare ~~~~~~~~~ //
        hp = ff.ButterworthHP(fs, fc, 0.7071);  // Prepare the high-pass filter
        lp = ff.ButterworthLP(fs, fs*0.25, 0.7071); // Prepare the low-pass filter
      }                                 // ~~~~~~~~~ Prepare ~~~~~~~~~ //
      inline T Process(T x) noexcept               // ~~~~~~~~~ Process ~~~~~~~~~ //
      {                                            // Process the input through both filters
        return lp.ProcessSample(hp.ProcessSample(x));
      }                                            // ~~~~~~~~~ Process ~~~~~~~~~ //
      inline void Clear(void) noexcept
      {                                    // ~~~~~~~~~ Clear ~~~~~~~~~ //
        hp.Clear();                        // Clear the high-pass filter
        lp.Clear();                        // Clear the low-pass filter
      }                                    // ~~~~~~~~~ Clear ~~~~~~~~~ //
  };
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Pluck pick bow position using our delay-branch
      // It features a integer delay, const Thiran fractional delay, variable Farrow fractional delay
      // to allow for vibrato and other modulation effects.
      // --------------------------------- //
    template<typename T=float,              // Our sample format
    size_t MaxLen=1<<15>                           // Thiran Filter order.
    class PositionEQ
    {
      size_t M{0};
      public:
        void Prepare(size_t looplen,double pos)
        {                           // ~~~~~~~~~ Prepare ~~~~~~~~~ //
            M=static_cast<size_t>(looplen*pos)&(MaxLen-1);
        }                           // ~~~~~~~~~ Prepare ~~~~~~~~~ //
      template<class Delay> inline T operator()(const Delay& dl) const noexcept
      {                              // ~~~~~~~~~ Peek operator ~~~~~~~~~ //
        return dl.Peek(0)-dl.Peek(M);// Get the difference between the current sample and the sample at position M
      }                              // ~~~~~~~~~ Peek operator ~~~~~~~~~ //
    };  
} // namespace sig::wg
