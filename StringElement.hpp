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
 * */
#pragma once
#include <cmath>
#include <random>
#include <cassert>
#include <experimental/simd>
#include <algorithm> // for std::min
#include "DelayBranch.hpp"
#include "DispersionAP.hpp"           
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
      double m=0.0,                      // Fractional delay in samples (Thiran)
      double m1=0.0,                     // Fractional delay in samples (Farrow)
      size_t o=K,                        // Order of the Farrow deinterpolator
      double p=0.15) noexcept            // Position in the string
    {                                    // ~~~~~~~~~ Prepare ~~~~~~~~~~ //
      if (fs <= 0.0 || lossfc>=0.5*fs) return false; // Sanitize sample rate
      if (f<=0.0||f>=22.5e3) return false;
      this->fs=fs;                       // Set the sample rate
      f0=f;                              // Set the fundamental frequency
      this->D=D;                         // Set the integer delay in samples
      this->mut=m;                        // Set the fractional delay in samples
      this->muf=m1;                      // Set the fractional delay for Farrow
      order=o;                           // Set the order of the Farrow deinterpolator
      pos=p;                             // Set the position in the string
      lfc=lossfc;                     // Set the loss filter cutoff frequency
      // Use provided integer delay for loop length
      size_t L=D?D:static_cast<size_t>(fs/(2.0*f0)+0.5f);
      // Initialize delay branches (ignore individual failures)
      fwd.Prepare(L,mut,muf);
      rev.Prepare(L,mut,muf);
      // Prepare the Delay Branch for the Group Delay
      const size_t gd=fwd.GroupDelay(L,o,o); // Compute the group delay for this structire
      for (size_t k=0;k<gd;++k)                 // For the length of the group delay....
      {                                 // Prime the Delay Branches
        fwd.Write(T(0));                // Write zero to the forward delay branch
        rev.Write(T(0));                // Write zero to the reverse delay branch
        fwd.Propagate(1);               // Propagate the zero sample through the forward delay branch
        rev.Propagate(1);               // Propagate the zero sample through the reverse delay branch
      }                                 // Done priming waveguide branches
      constexpr T eps=1e-5;             // Small epsilon value for numerical stability
      assert(std::fabs(fwd.Read())<eps&&std::fabs(rev.Read())<eps); // Ensure both branches are primed
      chain.Prepare(fs,alpha);           // Prepare the loss chain with sample rate and cutoff frequency
      stages=std::min<size_t>(st,disp.size());// Limit the number of dispersion stages to the size of the array
      for (size_t k=0;k<stages;++k)     // For each dispersion stage...
        disp[k].Prepare(alpha);          // Prepare each dispersion stage with the given coefficient
      bridge=ff.ButterworthHP(fs,fs*0.25,0.7071); // Prepare the bridge filter
      nut=ff.ButterworthHP(fs,40.0,0.7071); // Prepare the nut filter
      bridge=ff.ButterworthLP(fs,fs*0.45,0.7071); // Prepare the bridge filter
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
    // backward compatibility with higher-level code such as WaveguideVoice.
    inline bool Prep(double fs, double f, double lossfc,double q) noexcept
    {
      (void)q;     // Ignore the Q parameter for now, as it's not used in this context
      double lenlp=fs/(2.0*f);          // Compute looop length
      double id=std::floor(lenlp);      // Compute integer delay
      double fm=lenlp-id;                // Compute fractional delay
      mut=muf=fm; // Set fractional delay for Thiran and Farrow
      return Prepare(fs,f,lossfc,alpha,stages,id,mut,muf,order,pos);
    }
    void Excite(float s) noexcept
    {                                   // ~~~~~~~~~ Excite ~~~~~~~~~~ //
      fwd.Write(s);                    // Write the excitation sample to the forward delay branch
      fwd.Propagate(1);                // Propagate the sample through the forward delay branch
      rev.Write(-s);                   // Write the excitation sample to the reverse delay branch
      rev.Propagate(1);                // Propagate the sample through the reverse delay branch
    }                                   // ~~~~~~~~~ Excite ~~~~~~~~~~ //
    // Old API helpers.
    inline T Output(void) const noexcept
    {                                   // ~~~~~~~~~ Output ~~~~~~~~~~ //
      return static_cast<T>(pos);       // Return the current position in the string
    }                                   // ~~~~~~~~~ Output ~~~~~~~~~~ //
    void Propagate(size_t n) noexcept
    {
      for (size_t i=0;i<n;++i)          // For each sample to circulate.
      {                                 // Circulate ~~~~~~~~~~~~~~~~~~~~
        // Ensure branches have up-to-date fractional delay (vibrato, position mod)
        fwd.SetFractionalDelay(static_cast<T>(mut));
        fwd.SetMuFarrow(static_cast<T>(muf));
        rev.SetFractionalDelay(static_cast<T>(mut));
        rev.SetMuFarrow(static_cast<T>(muf));
        float pf=fwd.Read();                // Read through Thiran+Farrow
        float pr=rev.Read();                // Read through Thiran+Farrow
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
      prepared=false;                  // Set the prepared flag to false
      pos=0.0f;                        // Reset the position in the string
      f0=440.0f;                       // Reset the fundamental frequency
      D=0.0f;                          // Reset the integer delay in samples
      mut=0.0f;                         // Reset the fractional delay in samples
      muf=0.0f;                        // Reset the fractional delay for Farrow
      npwr=2e-4f;                      // Reset the noise power
      order=K;                         // Reset the order of the Farrow deinterpolator
    }                                   // ~~~~~~~~~ Clear ~~~~~~~~~~ //
    inline void SetPosition(double p) noexcept
    {
      pos=p;                            // Set the position in the string
    }
    inline double GetPosition(void) const noexcept
    {
      return pos;                       // Return the position in the string
      Prepare(fs,f0,lfc,alpha,stages,D,mut,muf,order,pos); // Re-prepare the string element with the new fundamental frequency

    }
    inline void SetPitchBendCents(double cents) noexcept
    {                                    // ~~~~~~~~~ SetPitchBendCents ~~~~~~~~~~ //
      double factor=std::pow(2.0,cents/1200.0);
      double newF0=f0*factor;
      if (newF0<=0.0) return;
      double loopLen=fs/(2.0*newF0);
      double id=std::floor(loopLen);
      double fm=loopLen-id;
      f0=newF0;
      D=id;
      mut=muf=fm;
      Prepare(fs,f0,lfc,alpha,stages,D,mut,muf,order,pos); // Re-prepare with new pitch bend
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
      Prepare(fs,f0,lfc,alpha,stages,D,mut,muf,order,pos); // Re-prepare the string element with the new fundamental frequency
    }
    inline void SetFundamental(double f) noexcept
    {
      f0=f;
      Prepare(fs,f0,lfc,alpha,stages,D,mut,muf,order,pos); // Re-prepare the string element with the new fundamental frequency
    }
    inline void SetIntegerDelay(double d) noexcept
    {
      D=d;
      Prepare(fs,f0,lfc,alpha,stages,D,mut,muf,order,pos); // Re-prepare with new integer delay
    }
    inline void SetFractionalDelay(double m) noexcept
    {
      mut=m;
      fwd.SetFractionalDelay(m);
      rev.SetFractionalDelay(m);
      Prepare(fs,f0,lfc,alpha,stages,D,mut,muf,order,pos); // Re-prepare with new fractional delay
    }
    inline void SetMuFarrow(double m1) noexcept
    {
      muf=m1;
      fwd.SetMuFarrow(m1);
      rev.SetMuFarrow(m1);
      Prepare(fs,f0,lfc,alpha,stages,D,mut,muf,order,pos); // Re-prepare the string element with the new fractional delay
    }
    inline void SetPositionEQ(double p) noexcept
    {
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
      Prepare(fs,f0,D,mut,muf,order,pos); // Re-prepare the string element with the new number of dispersion stages
    }
    inline size_t GetDispersionStages(void) const noexcept
    {
      return stages;                   // Return the number of dispersion stages
    }
    inline void SetDispersionCoefficient(double alpha) noexcept
    {
      for (size_t k=0;k<stages;++k)    // For each dispersion stage...
        disp[k].SetCoefficient(alpha); // Set the coefficient of each dispersion stage
      Prepare(fs,f0,D,mut,muf,order,pos); // Re-prepare the string element with the new number of dispersion stages
    }
    inline double GetDispersionCoefficient(void) const noexcept
    {
      return disp[0].GetCoefficient(); // Return the coefficient of the first dispersion stage
    }
    inline void SetBridgeFilter(double fc, double q) noexcept
    {
      bridge.SetCutoffFrequency(fc);    // Set the cutoff frequency of the bridge filter
      bridge.SetQualityFactor(q);       // Set the Q factor of the bridge filter
      Prepare(fs,f0,D,mut,muf,order,pos); // Re-prepare the string element with the new number of dispersion stages
    }
    inline void SetNutFilter(double fc, double q) noexcept
    {
      nut.SetCutoffFrequency(fc);       // Set the cutoff frequency of the nut filter
      nut.SetQualityFactor(q);          // Set the Q factor of the nut filter
      Prepare(fs,f0,D,mut,muf,order,pos); // Re-prepare the string element with the new number of dispersion stages
    }
    inline void SetOrder(size_t o) noexcept
    {
      order=o;
    }
    inline size_t GetOrder(void) const noexcept
    {
      return order;                    // Return the order of the Farrow deinterpolator
    }
    inline void SetLoopCutoff(double fc) noexcept
    {
        lfc=fc;                          // Set the loop cutoff frequency
        chain.Prepare(fs, lfc);       // Prepare the loss chain with the new cutoff frequency
        Prepare(fs,f0,D,mut,muf,order,pos); // Re-prepare the string element with the new number of dispersion stages
    }
    inline double GetLoopCutoff(void) const noexcept
    {
      return lfc;                       // Return the loop cutoff frequency
    }
    inline void SetLoopQ(double q) noexcept
    {
        lq=q;                             // Set the loop Q factor
        chain.Prepare(fs, lfc);       // Prepare the loss chain with the new cutoff frequency
        Prepare(fs,f0,D,mut,muf,order,pos); // Re-prepare the string element with the new number of dispersion stages
    }
    inline double GetLoopQ(void) const noexcept
    {
      return lq;                        // Return the loop Q factor
    }
    inline void SetAlpha(double a) noexcept
    {
      alpha=a;                          // Set the dispersion coefficient
      for (size_t k=0;k<stages;++k)    // For each dispersion stage...
        disp[k].SetCoefficient(alpha); // Set the coefficient of each dispersion stage
      Prepare(fs,f0,D,mut,muf,order,pos); // Re-prepare the string element with the new number of dispersion stages
    }
    inline double GetAlpha(void) const noexcept
    {
      return alpha;                    // Return the dispersion coefficient
    }
  private:
    Branch fwd{}, rev{};
    FilterFactory<T> ff; // Filter factory for creating filters
    LossChain<float> chain;             // Our loss chain
    std::array<DispersionAllPass<float>,2> disp{};
    size_t stages{2};                    // Number of dispersion stages
    BiQuad<float> bridge,nut;           // Bridge and nut filters
    PositionEQ<float> peq;              // Position EQ
    bool prepared{false};               // Preparation flag
    double lfc{6000.0};                 // Loss filter cutoff frequency
    double lq{0.7071};                 // Loss filter Q factor
    double fs{48000.0};                 // Sample rate
    double f0{440.0};                   // Fundamental frequency
    double pos{0.f};                    // Position in the string
    double alpha{0.18};                // Dispersion coefficient
    double D{0.0};                      // Integer delay in samples
    double mut{0.0};                    // Fractional delay in samples (Thiran)
    double muf{0.0};                    // Fractional delay in samples (Farrow)
    double npwr{2e-4f};                // Noise power
    size_t order{3};                    // Order of the Farrow deinterpolator
    std::minstd_rand rng; std::uniform_real_distribution<float> dist{-1.0f,1.0f}; // Random number generator for noise
    inline float Noise(void) noexcept
    {
      return dist(rng);                    // Return the generated noise
    }
    inline void StartNoise(void) noexcept
    {
      rng.seed(2025); // Seed the random number generator
    }    
};
}
