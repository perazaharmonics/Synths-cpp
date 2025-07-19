/* 
 * * Description:
 * * This file contains a delay line class used to represent a sample
 * * of sound travelling through a waveguide in either fwd (S+) or bwd (S-)
 * * direction. Two DelayBranch objects in opposite directions form one
 * * bidirectional tube/string. 
 * * 
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 * * SEQUENCE_SCENARIO:
 * *   The way to envision the plumbing here is...,
 * *   Imagine the following waveguide oversimplification:
 * *   You have a signal coming in from the left and need to
 * *   delay it by some samples, and then suddenly.... a fraction?!
 * *   the musician must be a virtuoso and is pitch bending...
 * *   We must propagate the left side incoming signal, to the right side
 * *   and then apply a fractional delay to it, so we can read it out
 * *   with a fractional delay back into a delay line encoded in size_t like
 * *   it was never REALLY hard in the first place.
 * *   So we apply the following sequence of filters:
 * * 
 * * SEQUENCE:
 * * 
 * * READ SIDE:
 * * WAVEGUIDE: ======================================================
 * * Signal -> DelayLine -> Thiran (mut) -> Farrow (muf[n]) -> Y
 * * WAVEGUIDE: ======================================================
 * * WRITE SIDE:======================================================
 * *  Signal <-DelayLine <- Thiran^-1 <- Farrow^-1 <- Y
 * * WAVEGUIDE: ======================================================
 * *
 * * 
 */
#pragma once
#include <algorithm>
#include <cstdio>
#include "DelayLine.hpp"
#include "FarrowFIR.hpp"
#include "ThiranAllPass.hpp"
#include "WGTypes.hpp"

#ifndef DBGP
#define DBGP(fmt, ...) std::printf("[DelayBranch] " fmt "\n", ##__VA_ARGS__)
#endif

namespace sig::wg
{
  template<typename T=float,              // Our sample format
    size_t MaxLen=1<<15,                  // Maximum length of the delay line in samples
    size_t K=3,                           // Farrow filter order
    size_t P=3>                           // Thiran Filter order.
  class DelayBranch final:public Node
  {
    public:
      DelayBranch(void) noexcept
      {
        dl=new sig::DelayLine<T,MaxLen>{}; // Create a new delay line with the specified maximum length.
        tip=new ThiranInterpolator<T,MaxLen,P>{}; // Create a new Thiran all-pass filter.
        tdip=new ThiranDeinterpolator<T,P>{}; // Create a new Thiran deinterpolator.
        fip=new FarrowInterpolator<T,MaxLen,K>{}; // Create a new Farrow interpolator.
        fdip=new FarrowDeinterpolator<T,MaxLen,K>{}; // Create a new Farrow deinterpolator.
      }
      ~DelayBranch(void) noexcept
      {
        delete dl;                      // Delete the delay line.
        dl=nullptr;
        delete tip;                     // Delete the Thiran all-pass filter.
        tip=nullptr;
        delete tdip;                    // Delete the Thiran deinterpolator.
        tdip=nullptr;
        delete fip;                     // Delete the Farrow interpolator.
        fip=nullptr;
        delete fdip;                    // Delete the Farrow deinterpolator.
        fdip=nullptr;
      }
    
    
    bool Prepare(
      size_t idelay,              // The integer delay for DelayLine.
      float m0,                         // The fractional delay for Thiran
      float m1=0.0) noexcept                // The fractional delay for Farrow
    {                                   // ----------- Prepare ------------- //
      if (idelay<1) return false;
      N=idelay;                         // Set integer delay in samples.
      if (N<K) return false;    /// Delay must be greater than filter order.
      mut=m0-std::floor(m0); // Set the fractional delay for Thiran.
      muf=m1-std::floor(m1); // Set the fractional delay for Farrow.
      dl->SetDelay(N+1);                // Set the delay line actual length;
      dl->Clear();                      // 
      tip->Prepare(N+mut, P);             // Prepare Thiran Interpolator
      tdip->Prepare(P, mut);            // Prepare Thiran Deinterpolator
      fip->SetOrder(K);                 // Set the order of the Farrow Interpolator
      fdip->SetOrder(K);                // Set the order of the Farrow Deinter
      fip->SetMu(muf);                // Set the fractional delay for Farrow Interpolator
      fdip->SetMu(-muf);               // Set the fractional delay for Farrow
      return true;                      // Return true if preparation was successful.
    }                                   // ----------- Prepare ------------- //
    // Runtime vibratio/pitch bend
    inline void SetMuFarrow(float m) noexcept
    {
      muf=m-std::floor(m); // Set the fractional delay for Farrow.
      fip->SetMu(muf);    // Set Thiran const fractional delay.
      fdip->SetMu(-muf); // Set the modulatable fractional delay for Farrow.
    }
    
    // Branch API matching StringElement
    inline size_t GetDelay() const noexcept { return N; }
    inline void SetFractionalDelay(float m) noexcept {
      mut = m - std::floor(m);
      tip->SetFractionalDelay(mut);
    }
    // ============================ Read ==============================
    // Read a sample from the delay line into fractional taps
    // This is the read side of the delay branch.
    // It reads the sample from the delay line, processes it through the Thiran all-pass
    // filter, and then processes it through the Farrow interpolator to get the output sample
    // ==================================================================
    //
    // READ SIDE: 
    // WAVEGUIDE: ======================================================
    // Signal -> DelayLine -> Thiran (mut) -> Farrow (muf[n]) -> Y
    // WAVEGUIDE: ======================================================
    // WRITE SIDE:
    // WAVEGUIDE: ======================================================
    //  Signal <-DelayLine <- Thiran^-1 <- Farrow^-1 <- Y
    // WAVEGUIDE: ======================================================
    //
    // =================================================================
    inline T Read(void)
    {
      // -------------------------------- //
      // 1. We are going backwards and reading from right to left.
      // Farrow needs all D...D-P taps in *original* form, so we first copy them
      // apply Thiran All Pass to each, then let Farrow do the mu1 mix.
      // -------------------------------- //
      // -------------------------------- //
      // Thiran stage: fill 'taps[]' 
      // -------------------------------- //
      DBGP("Read(): N=%zu  mut=%.6f  muf=%.6f", N, mut, muf);
      std::array<T,P+1> taps{};
      T x=T(0.f); // Initialize x to zero
      for (size_t k=0;k<=P;++k)           // For each Thiran coefficient
      {                                   // Circulate through Thiran's graph.
        // Check N-k_mut so it's always greater than zero
       double pos=static_cast<double>(k)+mut;
       T d=dl->ReadFrac(pos); // Read the fractional delay from the delay line.
       DBGP("  Thiran k=%zu pos=%.6f mu=%.6f  D=%.6f", k, pos, mut, d);
       x=tip->ProcessSample(d); // Process the sample through the Thiran all-pass filter.
       taps[k]=x;                        // Store the processed sample in the taps array.
       // K samples **behind** the head (negative offset)
       dl->WriteAt(static_cast<ptrdiff_t>(k),x);               // Write the processed sample back to the delay line.
       DBGP("  Thiran k=%zu  pos=%.3f  d=%.6f -> x=%.6f", k, pos, d, x);
      }                                   // Done setting Thiran taps.
      
      struct miniDL{
        const T* buf;
        T Peek(size_t i) const noexcept { return buf[i]; }
      } mini{taps.data()}; // Create a mini delay line with the Thiran taps.
      // -------------------------------- //
      // 2. Now we can process the taps through the Farrow interpolator.
      // -------------------------------- //
      T y{};                              // Output buffer
      T yF=fip->Process(dl,N/*D*/,&y);       // mu1[n] interpolation using actual delay line.
      T yT=x;                           // x is finished Thiran sample
      y+=T(0.5f)*(yT+yF);                 // The mixed signal.
      DBGP("[DelayBranch]   Thiran=%f  Farrow=%f  -> out=%f\n",
           double(yT), double(yF), double(y));
      return y;                           // Return the output sample.
    }
    // ============================ Write ==============================
    // This goes the other way around. The wave returns from the junction
    // so it actually gets sent back. 
    //
    // READ SIDE: 
    // WAVEGUIDE: ======================================================
    // Signal -> DelayLine -> Thiran (mut) -> Farrow (muf[n]) -> Y     
    // WAVEGUIDE: ======================================================
    // WRITE SIDE:                                                       
    // WAVEGUIDE: ======================================================
    //  Signal <-DelayLine <- Thiran^-1 <- Farrow^-1 <- Y
    // WAVEGUIDE: ======================================================
    //
    // This is the inverse of the Read function.
    // =================================================================
    void Write(T s) noexcept
    {
      DBGP("Write(): in=%.6f  N=%zu", s, N);
      // -------------------------------- //
      // 1. Inverse Farrow: distribute energy into delay line
      // ------------------------------- //
      DBGP("  pre-FDI  dl[N]=%.6f", dl->Peek(N));   // expect 0
      fdip->Process(s,*dl,N);            // Circulate through Farrow's graph.
      DBGP("  post-FDI dl[N]=%.6f", dl->Peek(N));
      // -------------------------------- //
      // 2. Inverse Thiran: phase-correct each tap
      // -------------------------------- //
      for (size_t k=0;k<=P;++k)           // For each tap in the Thiran filter
      {                                   // Circulate through Thiran's graph
        T x=dl->Peek(N-k);                // Get the sample at index N-k from the delay line.
        T y=tdip->ProcessSample(x);       // Process the sample through the Thiran deinterpolator.
        dl->WriteAt(static_cast<ptrdiff_t>(k),y);               // Write the processed sample back to the delay line.
        DBGP("  invThiran k=%zu   x=%.6f -> y=%.6f", k, x, y);
      }                                   // Done with the Thiran taps.
      haswritten=true;                    // We wrote something.
    }
    
    // Propagate the delay line by 'n' samples, normally n=1 per Tick()
    void Propagate(size_t n=1) noexcept
    {
    for (size_t i=0;i<n;++i)
    {
      const T before=dl->Read();
      dl->Write(before);
      DBGP("Propagate: hop %zu  sample=%.6f", i, before);
    }
    }
    void Clear(void) noexcept
    {
      dl->Clear();
    }                                   // Clear the delay line.
    inline size_t GetDelayLength(void) const noexcept
    {
      return N;               // Return the integer delay length in samples.
    }
    private:
      // ---------------------------------------------- //
      // The way to connects this is:
      // To read first write a signal to integer delay line, it'll be processed 
      // by the filter bank and returned through the Inverse DelayLine to be written out.
      // READ SIDE: 
      // WAVEGUIDE: ======================================================
      // Signal -> DelayLine -> Thiran (mut) -> Farrow (muf[n]) -> Y
      // WAVEGUIDE: ======================================================
      // WRITE SIDE:
      // WAVEGUIDE: ======================================================
      //  Signal <-DelayLine <- Thiran^-1 <- Farrow^-1 <- Y
      // WAVEGUIDE: ======================================================                                                                           <-

      // Members:
      // The base medium for the signal.
      DelayLine<float,MaxLen>* dl{nullptr};             // cache-friendly integer delay buffer
      bool haswritten{false};                            // flag for writes
      // Sub-sample shift for stretch-tuning and dispersion (Fixed mu = mut)
      ThiranInterpolator<T,MaxLen,P>* tip{nullptr};      // Thiran MF All-Pass interpolator.
      ThiranDeinterpolator<T,P>* tdip{nullptr}; // Thiran MF All-Pass deinterpolator. (phase corrector)
      // Sits after Thiran in the bank, handles very small fractional delays for vibrato swings
      // pitch bends with tiny polynomial error approximation. (up to 10 kHz)
      FarrowInterpolator<T,MaxLen,K>* fip{nullptr};       // Farrow Lagrange FIR interpolator
      FarrowDeinterpolator<T,MaxLen,K>* fdip{nullptr};  // Farrow Lagrange FIR deinterpolator ()
      // Integer and fractional parts of ideal delay line.
      /// Geometry related
      size_t N{0};
      float mut{0.0f}; // The fractional delay in samples for Thiran (fixed mu)
      
      // Modulation parameters
      float muf{0.0f}; // The fractional delay in sample for Farrow (variable mu) (modulation param)

      //static_assert(MaxLen > 0 && (MaxLen & (MaxLen - 1)) == 0, "MaxLen must be a power of two for efficient wrap-around.");
      //double len{MaxLen};             // The length of the delay line in samples.
   };
} // namespace sig::wg
