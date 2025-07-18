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
 */
#pragma once
#include <algorithm>
#include "DelayLine.hpp"
#include "FarrowFIR.hpp"
#include "ThiranAllPass.hpp"
#include "WGTypes.hpp"

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
        tip=new ThiranAllPass<T,MaxLen>{}; // Create a new Thiran all-pass filter.
        tdip=new ThiranDeinterpolator<T,MaxLen>{}; // Create a new Thiran deinterpolator.
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
      if (idelay<1||idelay>MaxLen-1) return false;
      N=idelay;                         // Set integer delay in samples.
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
      fip->SetMu(muf);
      fdip->SetMu(-muf); // Set the fractional delay for Farrow.
    }
    // Read: (DL -> Thiran -> Farrow)
    inline T Read(void)
    {
      // -------------------------------- //
      // 1. We are going backwards and reading from right to left.
      // Farrow needs all D...D-P taps in *original* form, so we first copy them
      // apply Thiran All Pass to each, then let Farrow do the mu1 mix.
      // -------------------------------- //
      std::array<T,P+1> taps{};           // Array to hold thiran taps.
      for (size_t i=0;i<=P;++k)           // For each tap in the Thiran filter
      {                                   // Circulate through Thiran's graph.
        T D=dl->Peek(N-k);                // Get the sample at index N-k from the delay line.
        taps[k]=tip->Write(D);            // Process the sample through the Thiran all-pass filter.
      }                                   // Done with the Thiran taps.
      // -------------------------------- //
      // 2. Now it gets strange because we need a special mini Delay-Line for Farrow
      // so we expose the tap array as a minimal ad-hoc "delay-line" to Farrow using
      // a lambda shim so that it gets trashed after we are done with the call.
      // -------------------------------- //
      struct mDL{                         // Fake delay line or intermediary
        const std::array<T,P+1>& buf;     // Member to hold the taps.
        T Peek(size_t i) const noexcept { return buf[idx]; }// Function to Peek
      } mini{taps};                       // Ad-hoc delay line for Farrow
      T y{};                              // Output buffer
      fip->Process(mini,0/*D*/,&y);       // mu1[n] interplation.
      return y;                           // Return the output sample.
    }
    // This goes the other way around. The wave returns from the junction
    // so it actually gets sent back. Farrow^-1 <- Thiran ^-1 <- DelayLine
    // This is the inverse of the Read function.
    void Write(T s) noexcept
    {
      // -------------------------------- //
      // 1. Inverse Farrow: Distributes energy across P+1 taps
      // So we first need to obtain provisional *integer* slot gains v[k] using
      // inverse Farrow
      // -------------------------------- //
      std::array<T,P+1> v{};              // Array to hold the provisional
      // Farrow taps.
      struct mDL{                         // Fake delay line or intermediary
        std::array<T,P+1>& buf;           // Member to hold the taps.
        void Write(size_t i, T x) noexcept { buf[i]=x; } // Function to write
        T Peek(size_t i) const noexcept { return buf[i]; } // Function to Peek
      } mini{v};                          // Ad-hoc delay line for Farrow
      fdip->Process(s,mini,0/*D*/);       // Process the sample through the Farrow deinterpolator.
      // -------------------------------- //
      // 2. Inverse Thiran: Now we pass it through inverse Thiran to correct phase shift
      // introduced by the Thiran Interpolator, then write it back to delay line so
      // we get our signal in the buffer fully interpolated.
      // -------------------------------- //
      for (int k=0;k<=P;++k)              // For each Thiran filter tap
        dl->WriteAt(N-k,tdip->Write(v[k])); // Write the sample to the delay line at index N-k.
      haswritten=true;                    // Mark that we have written a sample.
    }
    
    // Propagate the delay line by 'n' samples, normally n=1 per Tick()
    void Propagate(size_t n=1) noexcept
    {
      for (size_t i=0;i<n;++i)            // For each sample to propagate
        dl->Advance();                    // Hop along ring buffer.
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
      // The base medium for the signal.
      DelayLine<float,MaxLen>* dl{nullptr};             // cache-friendly integer delay buffer
      // ---------------------------------------------- //
      // The way to connects this is:
      // To read first write a signal to integer delay line, it'll be processed 
      // by the filter bank and returned through the Inverse DelayLine to be written out.
      // ----------------------------------------------- //
      //  Signal -> DelayLine -> Thiran (mut) -> Farrow (muf[n]) -> Signal Out
      //  Signal <- Farrow^-1 <- Thiran^-1 <- DelayLine <- Signal Out
      // ---------------------------------------------- //                                                                                <-
      // Sub-sample shift for stretch-tuning and dispersion (Fixed mu = mut)
      ThiranInterpolator<T,MaxLen,3>* th{nullptr};      // Thiran MF All-Pass interpolator.
      ThiranDeinterpolator<T,MaxLen,3>* thinv{nullptr}; // Thiran MF All-Pass deinterpolator. (phase corrector)
      // Sits after Thiran in the bank, handles very small fractional delays for vibrato swings
      // pitch bends with tiny polynomial error approximation. (up to 10 kHz)
      FarrowInterpolator<T,MaxLen,3>* f{nullptr};       // Farrow Lagrange FIR interpolator
      FarrowDeinterpolator<T,MaxLen,3>* finv{nullptr};  // Farrow Lagrange FIR deinterpolator ()
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
