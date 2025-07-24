 /*
 * *
 * * Filename: StringElement.hpp
 * *
 * * Description:
 * *  StringElement models an (ideal or damped) string element inside a
 * *  digital wave-guide network.
 * *  Two counter-propagating DelayBranches are coupled through a simple
 * *  one-pole loss filter.
 * * The string element can be excited with a single sample impulse
 * *  at its input.
 * * It has tunable parameters and modulation thanks to the Farrow deinterpolator,
 * * mulstiple dispersion stages, and a position EQ.
 * * Author:
 * *  JEP J. Enrique Peraza
 * *
  */
#pragma once
#include <cmath>
#include <random>
#include <experimental/simd>
#include <algorithm> // for std::min
#include "DelayBranch.hpp"
#include "DispersionAllPass.hpp"           
#include "FilterFactory.hpp"    
#include "WGFilters.hpp"
namespace sig::wg
{
//---------------------------------------------------------------------------
template<typename T=float,              // Our sample format
  size_t MaxLen=1<<15,                  // Maximum length of the delay line in samples
  size_t K=3,                           // Farrow filter order
  size_t P=3>                           // Thiran Filter order.                        // Thiran interpolator order
class StringElement final:public Node
{
  using Branch=DelayBranch<T,MaxLen,K,P>; // Our delay branch type
  public:
    bool Prepare(
      double fs,                        // Sample rate
      double f,                         // Fundamental frequency
      double lossfc=6000.0,             // Loss filter cutoff frequency
      double alpha=0.18,                 // Dispersion coefficient
      size_t st=2,                       // Number of dispersion stages
      double D=0.0,                      // Integer delay in samples
      double m=0.0,                      // Fractional delay in samples
      size_t o=K,                        // Order of the Farrow deinterpolator
      double p=0.15) noexcept            // Position in the string
    {                                    // ~~~~~~~~~ Prepare ~~~~~~~~~~ //
      if (fs <= 0.0 || lossfc>=0.5*fs) return false; // Sanitize sample rate
      this->fs=fs;                       // Set the sample rate
      f0=f;                              // Set the fundamental frequency
      this->D=D;                         // Set the integer delay in samples
      this->mu=m;                        // Set the fractional delay in samples
      order=o;                           // Set the order of the Farrow deinterpolator
      pos=p;                             // Set the position in the string
      // Use provided integer delay for loop length
      size_t L = static_cast<size_t>(D);
      // Initialize delay branches (ignore individual failures)
      fwd.Prepare(L, 0.0f, 0.0f);
      rev.Prepare(L, 0.0f, 0.0f);
      chain.Prepare(fs,alpha);           // Prepare the loss chain with sample rate and cutoff frequency
      stages=std::min<size_t>(st,disp.size());// Limit the number of dispersion stages to the size of the array
      for (size_t k=0;k<stages;++k)     // For each dispersion stage...
        disp[k].Prepare(alpha);          // Prepare each dispersion stage with the given coefficient
      bridge=ff.ButterworthHP(fs,fs*0.25,0.7071); // Prepare the bridge filter
      nut=ff.ButterworthHP(fs,40.0,0.7071); // Prepare the nut filter
      bridge=ff.ButterworthLP(fs,fs*0.45,0.7071); // Prepare the bridge filter
      deint.Prepare(D,mu,order);       // Prepare the Farrow deinterpolator with the given integer delay, fractional delay and order
      peq.Prepare(L,p);                 // Prepare the position EQ with the given loop length.
      StartNoise();                         // Generate noise for the string element
      prepared=true;                    // Set the prepared flag to true
      return true;                      // Return true if preparation was successful
    }                                   // ~~~~~~~~~ Prepare ~~~~~~~~~~ //
    // ---------------------------------------------------------------------
    // Extended preparation helpers
    //
    // Provide convenience overloads to compute internal delays automatically
    // or to call the extended Prepare with fewer arguments.  These ensure
    // backward compatibility with higher‑level code such as WaveguideVoice.
    inline bool Prep(double fs, double f, double lossfc,double q) noexcept
    {
      (void)q;     // Ignore the Q parameter for now, as it's not used in this context
      double lenlp=fs/(2.0*f);          // Compute looop length
      double id=std::floor(lenlp);      // Compute integer delay
      double fm=lenlp-id;                // Compute fractional delay
      return Prepare(fs, f, lossfc, 0.18, 2, id,fm,K,0.15);
    }
    void Excite(float s) noexcept
    {                                   // ~~~~~~~~~ Excite ~~~~~~~~~~ //
      fwd.Write(s);                    // Write the excitation sample to the forward delay branch
      fwd.Propagate(1);                // Propagate the sample through the forward delay branch
      rev.Write(-s);                   // Write the excitation sample to the reverse delay branch
      rev.Propagate(1);                // Propagate the sample through the reverse delay branch
    }                                   // ~~~~~~~~~ Excite ~~~~~~~~~~ //
    void Propagate(size_t n) noexcept
    {
      for (size_t i=0;i<n;++i)          // For each sample to circulate.
      {                                 // Circulate ~~~~~~~~~~~~~~~~~~~~
        float pf=deint.Process(fwd.Read()); // Process the forward delay branch sample through the Farrow deinterpolator
        float pr=deint.Process(rev.Read()); // Process the reverse delay branch sample through the Farrow deinterpolator
        float outb=chain.Process(pf);       // Process the forward sample through the loss filter
        for (size_t k=0;k<stages;++k)   // For each dispersion stage
          outb=disp[k].ProcessSample(-outb);// Process the output through the bridge filter
        float outn=chain.Process(pr);       // Process the reverse sample through the loss filter
        for (size_t k=0;k<stages;++k)   // For each dispersion stage
          outn=disp[k].ProcessSample(outn);// Process the output through the bridge filter
        outn=bridge.ProcessSample(-outn);// Process the output through the nut filter
        outb+=Noise()*npwr;            // Add noise to the output of the forward delay branch
        outn+=Noise()*npwr;            // Add noise to the output of the reverse delay branch
        fwd.Write(outb);               // Write the output of the forward delay branch
        fwd.Propagate(1);              // Propagate the output through the forward delay branch
        rev.Write(outn);               // Write the output of the reverse delay branch
        rev.Propagate(1);              // Propagate the output through the reverse delay branch
        pos=outb+outn;                 // Update the position in the string
      }                                // Circulate ~~~~~~~~~~~~~~~~~~~~
    }                                  // Propagate ~~~~~~~~~~~~~~~~~~~~
    void Clear(void) noexcept
    {                                   // ~~~~~~~~~ Clear ~~~~~~~~~~ //
      fwd.Clear();                     // Clear the forward delay branch
      rev.Clear();                     // Clear the reverse delay branch
      chain.Clear();                   // Clear the loss filter
      for (size_t k=0;k<stages;++k)    // For each dispersion stage...
        disp[k].Clear();               // Clear each dispersion stage
      bridge.Clear();                  // Clear the bridge filter
      nut.Clear();                     // Clear the nut filter
      deint.Clear();                   // Clear the Farrow deinterpolator
      prepared=false;                  // Set the prepared flag to false
    }                                   // ~~~~~~~~~ Clear ~~~~~~~~~~ //
    inline void SetPosition(double p) noexcept
    {
      pos=p;                            // Set the position in the string
    }
    inline double GetPosition(void) const noexcept
    {
      return pos;                       // Return the position in the string
    }
    // -----------------------------------------------------------------
    // Output helper
    //
    // Returns the current output of the string element.  The position
    // variable holds the sum of the counter‑propagating waves computed
    // during Propagate().  This method matches the Output() function
    // expected by higher‑level code such as WaveguideVoice.
    inline T Output(void) const noexcept
    {
      return static_cast<T>(pos);
    }
    // -----------------------------------------------------------------
    // Pitch bend helper
    //
    // Adjusts the fundamental frequency by the given number of cents
    // relative to the current f0.  A positive value bends the pitch up,
    // a negative value bends it down.  This recomputes the integer and
    // fractional delays and re‑prepares the element using the helper
    // Prepare overload that accepts (fs, f, D, mu, order, pos).
    inline void SetPitchBendCents(double cents) noexcept
    {                                    // ~~~~~~~~~ SetPitchBendCents ~~~~~~~~~ //
      // Convert cents to a frequency multiplier
      double factor=std::pow(2.0,cents/1200.0);
      double newF0=f0*factor;
      if (newF0<=0.0) return;
      // Compute new delays
      double loopLen=fs/(2.0*newF0);
      double id=std::floor(loopLen);
      double fm=loopLen-id;
      f0=newF0;
      D=id;
      mu=fm;
      Prepare(fs,f0,D,mu,order,pos);
    }
    inline double GetFundamental(void) const noexcept
    {
      return f0;                        // Return the fundamental frequency
    }
    inline double GetSampleRate(void) const noexcept
    {
      return fs;                        // Return the sample rate
    }
    inline void SetSampleRate(double s) noexcept
    {
      fs=s;
    }
    inline void SetFundamental(double f) noexcept
    {
      f0=f;
    }
    inline void SetIntegerDelay(double d) noexcept
    {
      D=d;
    }
    inline void SetFractionalDelay(double m) noexcept
    {
      mu=m;
    }
    inline void SetPositionEQ(double p) noexcept
    {
      // Re-prepare PositionEQ with current loop length
      size_t L = static_cast<size_t>(fs/(2.0*f0));
      peq.Prepare(L, p);
    }
    inline void SetNoisePower(double np) noexcept
    {
      npwr=np;                         // Set the noise power
    }
    inline double GetNoisePower(void) const noexcept
    {
      return npwr;                    // Return the noise power
    }
    inline void SetDispersionStages(size_t st) noexcept
    {
      stages=std::min<size_t>(st,disp.size()); // Set the number of dispersion stages
      Prepare(fs,f0,D,mu,order,pos);   // Re-prepare the string element with the new number of dispersion stages
    }
    inline size_t GetDispersionStages(void) const noexcept
    {
      return stages;                   // Return the number of dispersion stages
    }
    inline void SetDispersionCoefficient(double alpha) noexcept
    {
      for (size_t k=0;k<stages;++k)    // For each dispersion stage...
        disp[k].SetCoefficient(alpha); // Set the coefficient of each dispersion stage
      Prepare(fs,f0,D,mu,order,pos);   // Re-prepare the string element with the new dispersion coefficient
    }
    inline double GetDispersionCoefficient(void) const noexcept
    {
      return disp[0].GetCoefficient(); // Return the coefficient of the first dispersion stage
    }
    inline void SetBridgeFilter(double fc, double q) noexcept
    {
      bridge.SetCutoffFrequency(fc);    // Set the cutoff frequency of the bridge filter
      bridge.SetQualityFactor(q);       // Set the Q factor of the bridge filter
      Prepare(fs,f0,D,mu,order,pos);   // Re-prepare the string element with the new bridge filter parameters
    }
    inline void SetNutFilter(double fc, double q) noexcept
    {
      nut.SetCutoffFrequency(fc);       // Set the cutoff frequency of the nut filter
      nut.SetQualityFactor(q);          // Set the Q factor of the nut filter
      Prepare(fs,f0,D,mu,order,pos);   // Re-prepare the string element with the new nut filter parameters
    }
    inline void SetOrder(size_t o) noexcept
    {
      order=o;
    }
    inline size_t GetOrder(void) const noexcept
    {
      return order;                    // Return the order of the Farrow deinterpolator
    }
    inline void SetMu(double m) noexcept
    {
      mu=m;
    }
  private:
    Branch fwd{}, rev{};
    FilterFactory<T> ff; // Filter factory for creating filters
    // Filters:
    LossChain<float> chain;             // Our loss chain
    std::array<DispersionAllPass<float>,2> disp{};
    size_t stages{};                    // Number of dispersion stages
    BiQuad<float> bridge,nut;           // Bridge and nut filters
    FFIRDeinterpolator<float,MaxLen,K> deint; // Farrow deinterpolator
    PositionEQ<float> peq;              // Position EQ
    bool prepared{false};               // Preparation flag
    // State 
    double fs{48000.0};                 // Sample rate
    double f0{440.0};                   // Fundamental frequency
    double pos{0.f};                    // Position in the string
    double D{0.0};                      // Integer delay in samples
    double mu{0.0};                    // Fractional delay in samples
    double npwr{2e-4f};                // Noise power
    size_t order{3};                    // Order of the Farrow deinterpolator
    // Noise soirce
    std::minstd_rand rng; std::uniform_real_distribution<float> dist{-1.0f,1.0f}; // Random number generator for noise
    inline float Noise(void) noexcept
    {
      return dist(rng);                    // Return the generated noise
    }
    // Noise generation and Start
    inline void StartNoise(void) noexcept
    {
      rng.seed(2025); // Seed the random number generator
    }    
};
}