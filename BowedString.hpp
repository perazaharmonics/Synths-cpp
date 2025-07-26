  /*
 * * P&G Labs, LLC
 * *
 * * Filename: BowedString.hpp
 * * Description:
 * *   BowedString element (1-D wave-guide) driven by a simplified
 * *   Smith–Hausser elasto-plastic friction model.
 * *
 * *   Travelling-wave buffers are **DelayBranch** waveguides that observe
 * *   that can process a fractional delay encoded in the delay line:
 * *
 * *       read : DelayLine -> Thiran -> Farrow
 * *       write: Farrow⁻¹ <- Thiran⁻¹ <- DelayLine
 * *
 * * Author:
 * *   JEP  J. Enrique Peraza, P&G Labs, LLC
 * *
 */
#pragma once
#include <cassert>
#include <cmath>
#include <algorithm>

#include "WGTypes.hpp"
#include "DelayBranch.hpp"
#include "OnePole.hpp"
#include "WGFilters.hpp"
#include "DispersionAllPass.hpp"
#ifndef DBGP
#define DBGP(fmt, ...) std::printf("[BowedString] " fmt "\n", ##__VA_ARGS__)
#endif
namespace sig::wg
{
template<typename T=double,             // Our sample format
size_t MaxLen=1<<15,                    // Maximum length of the delay line in samples
size_t K=3,                             // Farrow filter order
size_t P=3>                             // Thiran Filter order.
class BowedString final:public Node
{
  // The most wonderful invention, thanks Julius O. Smith for your knowledge!
  // Reference: https://ccrma.stanford.edu/~jos/pasp/Delay_Branch_Waveguide_Networks.html
  using Branch=DelayBranch<T,MaxLen,K,P>;
public:
  // Output and settings API
  inline T Pressure(void) const noexcept { return static_cast<T>(out); }
  inline void SetLoopCutoff(double fc) noexcept { lfc = fc; chain.Prepare(fs, lfc); }
  inline void SetLoopQ(double q) noexcept { lq = q; }
  inline void SetBowPosition(double p) noexcept { pos = p; }
  
  bool Prepare(                         // Prepare the bowed string element.
    double fs,                          // Sample rate in Hz.
    double f,                           // The fundamental frequency in Hz.
    double pos=0.5,                     // The position of the bow on the string (0.0 to 1.0).
    double lossfc=6000.0,               // The cutoff frequency of the loss filter in Hz.
    double alpha=0.18,                  // The dispersion coefficient (0.0 to 0.5).
    size_t st=2,                        // Number of dispersion stages (1 or 2).
    double D=0.0,                       // The integer delay in samples (0.0 to fs/2).
    double m=0.0,                       // The fractional delay in samples (Thiran).
    double m1=0.0,                      // The fractional modulation delay in samples (Farrow).
    size_t o=K,                         // The order of the interpolation filters.
    double p=0.15) noexcept             // The position in the string (0.0 to 1.0).
  {                                     // ~~~~~~~~~~~~~~~ Prepare ~~~~~~~~~~~~~~~ //
    if(fs<=0.0||lossfc>=0.5*fs) return false;// Sanitize sample rate
    if(f<=0.0||f>=22050.0) return false;// Sanitize fundamental frequency
    if(D>0.0&&D<K) return false;        // This is a condition set forth by DSP theory.
    // Store user settings.             // ~~~~~~~~~~~~~~~~~~ //
    this->fs=fs;                        // Set the sample rate.
    f0=f;                               // Set the fundamental frequency.
    this->D=D;                          // Set the integer delay in samples.
    mut=m;                              // Set the fractional delay for Thiran.
    muf=m1;                             // Set the fractional modulation delay for Farrow.
    order=o;                            // Set the order of the interpolation filters.
    this->pos=p;                        // Set the position in the string.
    lfc=lossfc;                         // Set the cutoff frequency of the loss filter.
    this->alpha=alpha;                  // Set the dispersion coefficient.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Determine the length of the Delay Branch. We do this so we can
    // calculate the group delay of the delay line. This allows us to
    // prime the delay line such that the next time we read from it, the data
    // will be available in the first iteration of wave propagation
    // through the digital waveguide (collection of physically aware elements
    // and buffers that store the resulting wave energy).
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    size_t L=D?static_cast<size_t>(D):static_cast<size_t>(fs/(2.0*f)+0.5);
    fwd.Prepare(L,mut,muf);            // Prepare the forward delay branch.
    rev.Prepare(L,mut,muf);            // Prepare the reverse delay branch.
    // This is the priming loop that does what we described above.
    // It will prime the delay line with zeros. (We put zeroes a lot everywhere in DSP?)
    const size_t gd=fwd.GroupDelay(L,o,o);// The Group Delay of a Branch with these parameters.
    for(size_t i=0;i<gd;++i)          // For the "length" of the group delay....
    {                                 // Prime.
      fwd.Write(0.0);                 // Pump the fwd branch with zeroes
      fwd.Propagate(1);               // Propagate the wave.
      rev.Write(0.0);                 // Pump the rev branch with zeroes
      rev.Propagate(1);               // Propagate the wave.
    }                                 // Done adjusting for Group Delay.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Interpolation leads to approximation errors, which leads to 0s not being
    // entirely zero. See Valimaki's Doctoral's thesis for why is this so:
    // I suggest you read the whole thesis, it's amazing.
    // http://users.spa.aalto.fi/vpv/publications/vesan_vaitos/ch3_pt1_fir.pdf
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    constexpr T eps=1e-5;             // Tiny bit of energy we consider zero.
    //assert(std::fabs(fwd.Read())<eps&&std::fabs(rev.Read())<eps);
    DBGP("Branches have energy of: %.6f, %.6f", fwd.Read(), rev.Read());
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Each time a wave meets a junction, you get the standing wave phenomenon.
    // Part of the wave is reflected due to impedance mismatch, and another is
    // transmitted. The loss chain tries to simulate impedance mismatche loss
    // through various stages of filters. It is also why we have the following for
    // loop with the dispersion stages filters.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    chain.Prepare(fs, lfc);           // Prepare the freq dependent loss simulation filter
    bridge=ff.ButterworthHP(fs, fs*0.45, 0.7071); // Prepare the bridge junction filter
    nut=ff.ButterworthHP(fs, 40.0, 0.7071);    // Prepare the nut junction filter
    peq.Prepare(L, pos);                           // Prepare the position EQ filter
    prepared=true;                    // A loaded gun.
    return true;                      // Return true if preparation was successful.
  }                                   // ~~~~~~~~~~~~~~~ Prepare ~~~~~~~~~~~~~~~ //
  void ResetState(void) noexcept      // Clear the state machine.
  {                                   // ~~~~~~~~~~~~~~~ ResetState ~~~~~~~~~~~~~~~ //
    fwd.Clear();                      // Clear the forward delay branch.
    rev.Clear();                      // Clear the reverse delay branch.
    for(auto& d:disp)                 // For each dispersion stage...
      d.Clear();                      // Clear the dispersion all-pass filter.
    bridge.Clear();                   // Clear the bridge junction filter
    nut.Clear();                      // Clear the nut junction filter
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Clear loss chain and reflection filter
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    chain.Clear();                    // Clear the loss chain filters
    loss.Reset();                     // Reset the reflection filter
    vb=0.1f;                          // Reset the bow velocity to 0.1 m/s.
    pm=0.3f;                          // Reset the bow pressure to 0.3 N.
    pmax=10.f;                        // Reset the maximum bow pressure to 10 N.
    out=0.f;                          // Reset the output to 0.
    mut=0.0;                          // Reset the fractional delay for Thiran.
    muf=0.0;                          // Reset the fractional delay for Farrow.
    vph=0.0;                          // Reset the vibrato phase to 0.
    vd=0.0;                           // Reset the vibrato depth to 0.
    vf=0.0;                           // Reset the vibrato frequency to 0.
  }
  void Clear(void) noexcept
  {
    ResetState();                  // Clear the state machine.
  }
  // ~~~~~~~~~~~~~~~~~~~~ //
  // Propagate a signal through the waveguide n times.
  // This is the main processing function of the API.
  // ~~~~~~~~~~~~~~~~~~~~ //
  void Propagate(size_t n) noexcept
  {                                     // ~~~~~~~~~~~~~~~ Propagate ~~~~~~~~~~~~~~~ //
    for(size_t i=0;i<n;++i)             // For the number of samples to Tick() for..
    {                                   // Circulate wave through the waveguide model.
      UpdateMuLFO();                    // Modulate if we have to.
      double wf=fwd.Read();             // Read the forward wave.
      double wr=rev.Read();             // Read the reverse wave.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Process through the loss chain it could be air or whatever medium
      // the wave is travelling through.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      double lf=chain.Process(wf);      // Process the forward wave through the loss chain.
      double lr=chain.Process(wr);      // Process the reverse wave through the loss chain.
      for(size_t k=0;k<stages;++k)      // For each dispersion stage....
      {                                 // Pump through computational blocks.
        lf=disp[k].ProcessSample(lf);   // Process the forward wave through the dispersion all-pass filter.
        lr=disp[k].ProcessSample(lr);   // Process the reverse wave through the dispersion all-pass filter.
      }                                 // Done pumping through the dispersion stages.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Calculate the bow velocity and pressure waves
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      double vString=lf+lr;             // String velocity proxy
      double vRel=vString-vb;           // Bow-String relative velocity.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // The bow pressure is a function of the bow velocity and the bow pressure maximum.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // The bow pressure is the difference between the bow pressure and the string velocity.
      double dP=std::clamp(pm-vString,0.0,pmax); 
      double U=(dP>0.0&&dP<pmax)?S*mub*dP*std::pow(1.0-dP/pmax,2.0):0.0;
      double Pw=U*(rho*c/S);            // The bow pressure wave.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // The blend factor
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      double blend=0.5;                 // Blend factor for the bridge and nut junctions.
      double ofw=(1.0-blend)*lf-blend*wf; // Forward wave at the bridge junction
      double ore=(1.0-blend)*lr-blend*wr; // Reverse wave at the bridge junction
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Now write the waves to the delay branches.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      fwd.Write(ofw);                   // Write the forward wave to the forward delay branch.
      fwd.Propagate(1);                 // Propagate the wave.
      rev.Write(ore);                  // Write the reverse wave to the reverse delay branch.
      rev.Propagate(1);                 // Propagate the wave.
    }                                   // Done propagating the wave through the waveguide.
  }                                     // ~~~~~~~~~~~~~~~ Propagate ~~~~~~~~~~~~~~~ //
  // Bow an impulse to the string model.
  void Bow(T v) noexcept                // 
  {                                     // ~~~~~~~~~~~~~~~ Bow ~~~~~~~~~~~~~~~ //
    fwd.Write(v);                       // Write the bow impulse to the forward delay branch.
    rev.Write(-v);                      // Write the bow impulse to the reverse delay branch.
  }                                     // ~~~~~~~~~~~~~~~ Bow ~~~~~~~~~~~~~~~ //
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  // Modulation API:
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  // Enable and disable modulation.
  inline void EnableVibrato(double depth,double rate){vd=depth;vf=rate;}
  inline void DisableVibrato(void){vd=0.0;}
  // Set Thiran Fractional delay
  inline void SetFractionalDelay(double m) 
  {                                     // ~~~~~~~~~~~~~~~ SetFractionalDelay ~~~~~~~~~~~~~~~ //
    mut=m;                              // Set the fractional delay for Thiran.
    fwd.SetFractionalDelay(m);          // Set the fractional delay for Thiran.
    rev.SetFractionalDelay(m);          // Set the fractional delay for Thiran.
    // Prepare the branches with the new fractional delay.
    Prepare(fs,f0,pos,lfc,alpha,stages,D,mut,muf,order,pos);
  }
  inline void SetIntegerDelay(double d) // Set the integer length of the delay line.
  {                                     // 
    D=d;                                // Set the integer delay in samples.
    Prepare(fs,f0,pos,lfc,alpha,stages,D,mut,muf,order,pos);
  }
  inline void SetMuFarrow(double m1)    // Set fracional delay for modulation Farrow FIR
  {                                     // ~~~~~~~~~~~~~~~ SetMuFarrow ~~~~~~~~~~~~~~~ //
    muf=m1;                             // Set the fractional delay for modulation Farrow.
    fwd.SetMuFarrow(m1);                // Set the fractional delay for modulation Farrow.
    rev.SetMuFarrow(m1);                // Set the fractional delay for modulation Farrow.
    // Prepare the branches with the new fractional delay.
    Prepare(fs,f0,pos,lfc,alpha,stages,D,mut,muf,order,pos);
  }
  inline void SetOrder(size_t o){order=o;}
  // Enabled thanks to mu farrow.
  inline void PitchBend(double cents)   // Pitch bend the string by a number of cents.
  {                                     // ~~~~~~~~~~~~~~~ PitchBend ~~~~~~~~~~~~~~~ //
    double fac=std::pow(2.0,cents/1200.0);// Calculate the pitch bend factor.
    Prepare(fs,f0*fac,pos,lfc,alpha,stages,D,mut,muf,order,pos);
  }                                     // ~~~~~~~~~~~~~~~ PitchBend ~~~~~~~~~~~~~~~ //
private:
  void UpdateMuLFO(void) noexcept       // Modulate according to mu farrow.
  {                                     // ~~~~~~~~~~~~~~~ UpdateMuLFO ~~~~~~~~~~~~~~~ //
    if(vd<=0.0||vf<=0.0||vf>=0.5*fs) return;// If no mod parameters do nothing....
    vph=std::fmod(vph*2*M_PI*vf/fs,2*M_PI);// Calculate the vibrato phase.
    double d=vd*std::sin(vph);          // Calculate the vibrato depth.
    double mt=std::clamp(mut+d,-0.95,0.95);// Recalculate Thiran fractional delay [-0.95,0.95]
    double mf=std::clamp(muf+d,-0.95,0.95);// Recalculate Farrow fractional delay [-0.95,0.95]
    fwd.SetFractionalDelay(mt);        // Set the fractional delay for Thiran.
    rev.SetFractionalDelay(mf);        // Set the fractional delay for Thiran.
    fwd.SetMuFarrow(mf);               // Set the fractional delay for modulation Farrow.
    rev.SetMuFarrow(mf);               // Set the fractional delay for modulation Farrow.
  }                                    // ~~~~~~~~~~~~~~~ UpdateMuLFO ~~~~~~~~~~~~~~~ //
private:
  // Forward and backward propagation waveguides.
  Branch fwd{},rev{};                   // Waveguide reflection branches.
  OnePole<T>loss;                       // Reflection filter for the waveguide.
  LossChain<T>chain;                    // Loss chain for the waveguide.
  std::array<DispersionAllPass<T>,2>disp{};// Dispersion stages for the waveguide.
  BiQuad<T>bridge,nut;                  //Nut junction and bridge junction filters.
  FilterFactory<T>ff;
  PositionEQ<T>peq;

  double fs{48000.0},f0{440.0},lfc{5000.0},lq{0.707},pos{0.1},vb{0.1},pm{0.3},
         pmax{10.0},S{1e-4},mub{1.0},rho{1.21},c{343.0},mut{0.0},muf{0.0},vph{0.0},
         vd{0.0},vf{0.0},out{0.0},alpha{0.18};

  size_t D{0},order{K},stages{2};
  bool prepared{false};
};
}
