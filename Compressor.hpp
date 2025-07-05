/*
 * 
 *
 * Filename: Compressor.hpp
 * Description:
 * Wrapper for the FaustCompressorDsp Faust-based Compressor effect.
 * Applies the same Compressor separately to left and right channels.
 *
 * Author:
 * JEP J. Enrique Peraza
 */
#pragma once
#include "AudioBuffer.hpp"
#include "FilterFactory.hpp"
#include <cmath>
namespace sig::effects
{
  template<typename T = float>
  class Compressor
  {
    public:
    Compressor(void) noexcept=default;
    ~Compressor(void) noexcept=default;
    // Must be called before processing, and the API tuning parameters
    // must be called before calling Prepare().
    bool Prepare(double fs, size_t /*blocksiz*/) noexcept
    {
      this->fs=fs;
      static bool envInitialized=false; // Flag to check if the envelope detector is initialized.
      if (!envInitialized) // If the envelope detector is not initialized...
      {
        envDetector=filterFactory.EnvelopeAttackDetector(fs,attackMs); // Create an envelope detector for attack.
        envInitialized=true; // Set the flag to true to indicate that the envelope detector is initialized.
      }
      static bool releaseInit=false;
      if (!releaseInit) // If the release detector is not initialized...
      {
        releaseDetector=filterFactory.EnvelopeReleaseDetector(fs,releaseMs); // Create an envelope detector for release.
        releaseInit=true; // Set the flag to true to indicate that the release detector is initialized.
      }
      return true;
    }
    bool Process(AudioBuffer<T>& io, size_t n) noexcept
    {
      auto* L=io.Channel(0);
      auto* R=io.Channel(1);
      if (!L||!R) return false;
      const T linThresh=std::pow(10.f,thresholdDb/20.f);
      for (size_t i=0;i<n;++i)
      {
        // Side-chan level (peak detector)
        T inL=std::fabs(L[i]),inR=std::fabs(R[i]);
        T level=std::max(inL,inR);
        // Envelope detector with attach and release
        if (level>env)
          env=envDetector.ProcessSample(level);
        else
          releaseDetector.ProcessSample(level);
        // Compute gain factor
        T gain=1;
        if (env>linThresh)
        {
          T dbIn=20*std::log10(env);
          T dbOut=thresholdDb+(dbIn-thresholdDb)/ratio;
          gain=std::pow(10.f,(dbOut-dbIn)/20.f);
        }
        // Mix calculations into buffers
        L[i]*=gain;
        R[i]*=gain; // Apply gain to both channels.
      }
      return true;
    }
    // API
    void SetThresholdDB(T dB) noexcept { this->thresholdDb=dB; }
    void SetRatio(T r) noexcept {this->ratio=r; }
    void SetAttackMs(T ms) noexcept { attackMs=ms; envDetector=filterFactory.EnvelopeAttackDetector(fs,ms); }
    void SetReleaseMs(T ms) noexcept { releaseMs=ms; releaseDetector=filterFactory.EnvelopeReleaseDetector(fs,ms); }
    T GetThresholdDB(void) const noexcept { return thresholdDb; }
    T GetRatio(void) const noexcept { return ratio; }
    T GetAttackMs(void) const noexcept { return attackMs; }
    T GetReleaseMs(void) const noexcept { return releaseMs; }
    inline double GetSampleRate(void) const noexcept { return fs; } // Get the sample rate for the compressor effect.
    inline void SetSampleRate(double fs) noexcept { this->fs=fs; } // Set the sample rate for the compressor effect.

    // Reset the internal state of the compressor effect.
    void Reset(void) noexcept
    {
      env=0; // Reset the envelope detector state.
      envDetector.Reset(); // Reset the envelope detector.
      releaseDetector.Reset(); // Reset the release detector.
    }
    private:
      double fs=48000.0; // Sample rate.
      T thresholdDb=-24.0f; // Threshold in dB.
      T ratio=4.0f; // Compression ratio.
      T attackMs=10.0f; // Attack time in milliseconds.
      T releaseMs=100.0f; // Release time in milliseconds.
      T env=0.f; // Current envelope level.
      OnePole<T> envDetector; // Envelope detector for attack.
      OnePole<T> releaseDetector; // Envelope detector for release.
      FilterFactory<T> filterFactory;
  };
}
