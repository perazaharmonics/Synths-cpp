/*
* * P&G Labs
* * 
* * Filename: FaustRingModulationP.hpp
* * Description:
* * This is actually some really old code that I wrote in Faust using CCRMA's
* * Computer Music lectures, and Julius O. Smith's book collection. It has
* * been exported to C++ using the Faust export tool. Then I just fed it to
* * a machine to repurpose it for the PGSynth project.
* *
* * Author:
* * JEP J. Enrique Peraza
* *
*/
#pragma once

#include "AudioBuffer.hpp"
#include <vector>

namespace sig::effects {

  template<typename T=float>
  class RingModulation {
  public:
    RingModulation(void) noexcept = default;
    ~RingModulation(void) noexcept = default;
    RingModulation(const RingModulation&) = delete;
    RingModulation& operator=(const RingModulation&) = delete;

    // Prepare engine with sample rate and block size
    bool Prepare(double sampleRate, size_t blockSize) noexcept 
    {
      this->fs=sampleRate;
      this->bs=blockSize;
      this->phaseincr=float(2*M_PI*modfr/fs);
      return true;
    }

    // Process interleaved or single-channel buffer through ring modulation + echo
    bool Process(AudioBuffer<T>& io, size_t n) noexcept 
    {
      auto* L=io.Channel(0); // Pointer to the left channel buffer.
      auto* R=io.Channel(1); // Pointer to the right channel buffer.
      if (L==nullptr || R==nullptr) return false; // Check if the buffers
      for (size_t i=0;i<n;++i) // For each sample in the
      {                                 // Ring modulate...
        // LFO
        float osc=std::sin(phase);      // Our modulation oscillator.
        phase+=phaseincr;               // Increment the phase.
        if (phase>=2*M_PI) phase-=2*M_PI; // Wrap the phase
        // Feedback ring-mod: mix previous out back into input
        T inL=L[i]+prevL*feedb;
        T inR=R[i]+prevR*feedb; // Mix previous output with current input.
        // Apply ring modulation
        T modL=inL*(1.0f-depth+depth*osc);
        T modR=inR*(1.0f-depth+depth*osc); // Apply modulation depth.
        // Store for next feedback iteration
        prevL=modL;
        prevR=modR;
        // Mix dry + wet into the buffers.
        L[i]=dry*L[i]+wet*modL; 
        R[i]=dry*R[i]+wet*modR; 
      }                                 // End of sample processing loop.
      return true;                      // Return true if processing is successful.
    }
    // Reset the internal state of the ring modulation effect.
    void Reset(void) noexcept 
    {
      phase = 0; // Reset the phase accumulator.
      prevL = T(0); // Reset previous left sample.
      prevR = T(0); // Reset previous right sample.
    }
    // API: Getters and Setters, always call before Prepare()
    inline double GetSampleRate(void) const noexcept { return fs; } // Get the sample rate for the ring modulation effect.
    inline void SetSampleRate(double fs) noexcept { this->fs=fs; } // Set the sample rate for the ring modulation effect.
    inline size_t GetBlockSize(void) const noexcept { return bs; } // Get the block size for the ring modulation effect.
    inline void SetBlockSize(size_t bs) noexcept { this->bs=bs; } // Set the block size for the ring modulation effect.
    inline void SetChannels(size_t nch) noexcept { this->nch=nch; } // Set the number of channels for the ring modulation effect.
    inline size_t GetChannels(void) const noexcept { return nch; } //

    // Ring-Mod API:
    inline void SetModulationFreq(float freq) noexcept { modfr=freq; phaseincr=float(2*M_PI*modfr/fs); } // Set the modulation frequency for the ring modulation effect.
    inline float GetModulationFreq(void) const noexcept { return modfr; } // Get the modulation frequency for the ring modulation effect.
    inline void SetDepth(float d) noexcept { depth=d; } // Set the modulation depth for the ring modulation effect.
    inline float GetDepth(void) const noexcept { return depth; } // Get the modulation depth for the ring modulation effect.
    inline void SetFeedback(float fb) noexcept { feedb=fb; } // Set the feedback amount for the ring modulation effect.
    inline float GetFeedback(void) const noexcept { return feedb; } // Get the feedback amount for the ring modulation effect.
    inline void SetDry(float d) noexcept { dry=d; } // Set the dry level for the ring modulation effect.
    inline float GetDry(void) const noexcept { return dry; } // Get the dry level for the ring modulation effect.
    inline void SetWet(float w) noexcept { wet=w; } // Set the wet level for the ring modulation effect.
    inline float GetWet(void) const noexcept { return wet; } // Get the wet level for the ring modulation effect.
    inline void SetPhase(float p) noexcept { phase=p; } // Set the phase for the ring modulation effect.
    inline float GetPhase(void) const noexcept { return phase; } // Get the current phase
    inline float GetPhaseIncrement(void) const noexcept { return phaseincr; } // Get the phase increment for the ring modulation effect.
  private:
    double fs{48000.0};
    size_t bs{512}, nch{2};
    

    // Ring mod parameters
    float modfr{800.0f};    // Modulation frequency in Hz
    float depth{1.0f};       // Modulation depth
    float feedb=0.0f;        // Feedback amount
    float dry=0.0f;          // Dry level
    float wet=1.0f;          // Wet level
    float phase=0;           // Phase accumulator for modulation
    float phaseincr=0;       // Phase increment for modulation frequency
    T prevL=T(0), prevR=T(0); // Previous samples for left and right channels
  };

} // namespace sig::effects
