 /*
 * *
 * * Filename: SDFDN.hpp
 * * 
 * * Description:
 * *  A Spectral Diffusion Feedback Delay Network (SDFDN) is a specialized application of a Feedback Delay Network
 * *  (FDN) design to create spectral diffusion within the mesh of the FDN. Branches of the FDN are modulated
 * *  using a VibLFO, which is a low-frequency oscillator (LFO) that modulates the delay lines in the FDN.
 * *  This VibLFO can also be wave-shaped to turn the LFO into a noise source, and feed noise into the FDN,
 * *  or layers of the FDN connected in series to a DelayBranch waveguide hosting the carrier which is about to
 * *  enter into a noisy, maybe high reflection environment. 
 * *
 * *  Author:
 * *    JEP  J. Enrique Peraza
 */
#pragma once
#include<cmath>
#include<vector>
#include<array>
#include<random>
#include<algorithm>
#include<atomic>
#include<cstddef>
#include"FeedbackDelayNetwork.hpp"
#include"VibLFO.hpp"
#include"Waveshaper.hpp"

namespace sig::wg{
/*---------------------------- Helpers ---------------------------*/
namespace detail{template<typename T>constexpr T TWOPI(){return static_cast<T>(6.2831853071795864769);} }

/*-------------------- Spectral Diffuser FDN ---------------------*/
template<typename T=float,                // Data processing type.
  size_t MaxLen=1<<15,                    // Maxlen of the delay line in samples.
  size_t K=3,                             // Farrow filter order
  size_t P=3,                             // Thiran filter order
  size_t N=8>                             // Number of taps in the FDN.
class SDFDN 
{
    static_assert(N>=4,"SpectralDiffuserFDN requires =4 taps");
    using FDN=FeedbackDelayNetwork<T,MaxLen,K,P,N>;
public:
    SDFDN (void)noexcept:rng(0xCEFAFD)
    {                                   // ~~~~~~~~~~ Constructor ~~~~~~~~~~ //
      constexpr T phi=.6180339887498948;// The golden ratio.
      for(size_t i=0;i<N;++i)           // For the number of taps in the FDN.
      {                                 // Initialize the LFOs.
        lfo[i].Enable(ld,lr);           // Enable the LFO with depth and rate
        lfo[i].SetPhase(std::fmod(phi*i,T(1))*detail::TWOPI<T>());
        dg[i]=1;                        // Set the depth gain for each tap
      }                                 // Initialize the FDN.
    }                                   // ~~~~~~~~~~ Constructor ~~~~~~~~~~ //
    ~SDFDN(void)noexcept = default;     // ~~~~~~~~~~ Destructor ~~~~~~~~~~ //
    bool Prepare (double fs,size_t bs) noexcept
    {                                   // ~~~~~~~~~~ Prepare ~~~~~~~~~~~~~~~ //
      // Sanitize input! Prevent over Nyquist limit Aliasing!!!!!!!!  
      if(fs<=0||fs>=0.5*fs||bs<1)return false;     
      this->fs=fs;                      // Set the sample rate
      this->bs=bs;                      // Set the block size
      std::array<double,N>d{};          // Default delays in seconds
      for(size_t i=0;i<N;++i)           // For each tap
        d[i]=.012+.004*std::pow(double(i),1.18);// Set the default delay
      fdn.SetDelays(d);                 // Set the default delays in the FDN
      fdn.SetFeedbackMatrix(detail::MatrixKind::Hadamard);// Set the feedback matrix to Hadamard
      fdn.SetWetMix(1);                 // Set the wet mix to 1.0
      std::array<double,N>lp{},sh{};    // Lowpass and shelf cutoffs
      for(size_t i=0;i<N;++i)           // For each tap
      {                                 // Set the lowpass and shelf -3db cutoff frequencies
       lp[i]=fs*0.45;                   // Lowpass cutoff at 45% of the sample rate
       sh[i]=7000+300*i;                // Shelf cutoff at 7kHz + 300Hz per tap
      }                                 // Done setting up FDN filters.
      bshelf=sh;                        // Remember we are using these shelf cutoffs.
      fdn.SetDamperCutoffs(lp);         // Set the damper cutoffs for the FDN
      fdn.SetShelfCutoffs(sh);          // Set the shelf cutoffs for the FDN
      fdn.Prepare(fs,bs);               // Prepare the FDN with the sample rate and block size
      return true;                      // All good if we got here.
    }                                   // ~~~~~~~~~~ Prepare ~~~~~~~~~~~~~~~ //
    inline bool Process (                // Process a slice of frames from the block.
      const T* in,                      // Pointer to input stream.
      T* const oL,                      // Where to store the processed signal left channel.
      T* const oR,                      // Where to store the processed signal right channel.
      size_t nFrames) noexcept          // The number of frames per block to process.
    {                                   // ~~~~~~~~~~ Process ~~~~~~~~~~~~~~~ //
      Tick();                           // Update the FDN and LFOs to next state.
      return fdn.Process(in,oL,oR,nFrames);// Forward the signal to the FDN and get output.
    }                                   // ~~~~~~~~~~ Process ~~~~~~~~~~~~~~~ //
    inline void Clear(void) noexcept
    {                                   // ~~~~~~~~~~~ Clear ~~~~~~~~~~~~~~~ //
      fdn.Clear();                      // Clear the FDN state
      for(auto&l:lfo)                   // For every LFO
        l.Clear();                      // Clear the LFO state
    }                                   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // API:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    inline void SetLFOShapeID (int id) noexcept 
    {                                   // ~~~~~~~~~~ SetLFOShapeID ~~~~~~~~~~ //
        sid.store(id%20);               // We want this waveform.
    }                                   // ~~~~~~~~~~ SetLFOShapeID ~~~~~~~~~~ //
    inline void SetTapDepths (const std::array<T,N>&d) noexcept 
    {
        dg=d;                           // Set the tap depths gain.
    }
    inline void SetDetuneAmount (T d) noexcept 
    {
        dtn=Clamp(d);                   // Set the detune amount
    }             
    inline void SetLFORate (T f) noexcept
    {
        f=f<.01?.01:f;                  // Clamp the LFO rate to a minimum of 0.01Hz
        for(auto&l:lfo)                 // For each LFO
            l.SetRate(f);               // Set the LFO rate
    } 
    inline void SetShimmerAmount (T v) noexcept 
    {
        shmr=Clamp(v);                  // Set how much shimmer we want
    }           
    inline void SetCloudyAmount (T v) noexcept 
    {
        cld=Clamp(v);                   // Set how much cloudiness we want
    }            
    inline void SetNoiseDensity (T v) noexcept 
    {
        dns=Clamp(v);                   // Set how much noise density we want
    }            
    inline void Freeze (bool f) noexcept 
    {
        freeze.store(f);                // If to freeze the FDN or not.
    }               
    inline bool IsFrozen (void) const noexcept 
    {
        return freeze.load();          // True if we froze the FDN
    }    
    inline void SetWetDryMix (T mix) noexcept 
    {
        fdn.SetWetMix(Clamp(mix));      // Set the wet/dry mix for the FDN
    }
    inline T GetWetDryMix (void) const noexcept 
    {
        return fdn.GetWetMix();         // Get the wet/dry mix from the FDN
    } 
    inline void SetSampleRate (double s) noexcept 
    {
        fs=s;
        fdn.SetSampleRate(s);
        for(auto&l:lfo)
          l.SetSampleRate(s);
    }
    inline double GetSampleRate(void) const noexcept {return fs;}
    inline void SetPreDelay(double s) noexcept
    {
      predelay.Prepare(static_cast<size_t>(s),0.0f,0.0f); 
    }
    void SetFractionalDelay(const std::array<double,N>& s) noexcept
    {                                   // Set Thiran fractional delays
      fdn.SetFractionalDelay(s);        // Set the Thiran fractional delays in the FDN
    }                                   // Set Thiran fractional delays
    void SetMuFarrow(const std::array<double,N>& s) noexcept
    {                                   // Set Farrow fractional delays.
      fdn.SetMuFarrow(s);               // Set the Farrow fractional delays in the FDN
    }                                   // Set Farrow fractional delays.
    inline void SetDelays(const std::array<double,N>& s) noexcept
    {
      fdn.SetDelays(s);
    }
private:
    // Utilities
    // Clamp a value between 0 and 1
    static constexpr T Clamp(T x) noexcept { return x<T(0)?T(0):x>T(1)?T(1):x; }
    // Block-wise update
    void Tick (void) noexcept
    {                                   // ~~~~~~~~~~ Tick ~~~~~~~~~~~~~~~ //
      std::array<double,N>mf{};         // Farrow modulation fractional delay
      // Our Uniform distribution (rectangle)
      T c=cld.load();std::uniform_real_distribution<T>dist(-1,1);
      for(size_t i=0;i<N;++i)           // For each tap.
      {                                 // Update the LFO and wave-shape it.
        T s=lfo[i].Tick(fs);            // Get the LFO value
        // Wave-shape the LFO value     //
        T v=WaveShaper<T>::Shape(s,lfo[i].GetPhase(),sid.load(),rng)*dg[i];
        v+=dtn.load()*dist(rng)*.5;     // Add detune to the LFO value
        mf[i]=.5*double(c)*v;           // Scale the LFO value by the cloudiness and depth gain
      }                                 // Done waveshaping the LFO.
      fdn.SetMuFarrow(mf);              // Set the Farrow modulation fractional delays in the FDN
      auto m=fdn.GetFeedbackMatrix();   // What feedback matrix are we using? Get that.
      T d=dns.load();                   // Load the noise density desired.
      for(auto&row:m)                   // For each row in the feedback matrix
        for(auto&v:row)                 // For each column vector in the row
          v*=d;                         // Scale the feedback matrix by the noise density
      fdn.SetFeedbackMatrix(m);         // Set the modified feedback matrix in the FDN
      std::array<double,N>cf{};         // Cutoff frequencies for the low-shelf filters
      T sm=shmr.load();                 // Shimmer amount
      for(size_t i=0;i<N;++i)           // For each tap
        cf[i]=bshelf[i]+6000*sm;        // Set the low-shelf cutoff frequency
      fdn.SetShelfCutoffs(cf);          // Set the low-shelf cutoff frequencies in the FDN
      if(freeze.load())                 // Were our parameters frozen?
      {                                 // Yes, so set feedback matrix to identity.
        auto I=detail::identity<T,N>(); // Create an identity matrix
        fdn.SetFeedbackMatrix(I);       // Set the feedback matrix to identity
      }                                 // Now everything should be as it was before freezing.
    }                                   // ~~~~~~~~~~~~~~ Tick ~~~~~~~~~~~~~~~~~ //
private:
    // Members
    FDN fdn;                           // Our base class, the FDN we operate on.
    std::array<VibLFO<T>,N>lfo;        // LFOs to modulate the delay lines
    std::array<T,N>dg{};               // Depth gain for each tap
    std::array<double,N>bshelf{};      // Low-shelf cutoff frequencies for each tap
    std::atomic<int>sid{0};            // Waveshaper wave ID.
    std::atomic<T>shmr{.5},cld{.4},dns{.9},dtn{.1};// Shimmer, cloudiness, noise density, detune parameters
    std::atomic<bool>freeze{false};    // True if user froze the FDN parameters.
    double fs{48000.0};                  // Sample rate
    double lr{.15};                    // LFO rate
    double ld{.5};                     // LFO depth
    double lp{0.0};                    // LFO phase
    size_t bs{256};                    // Our block processing size.
    std::mt19937 rng;                  // Our random number generator
};

}/*namespace sig::wg*/