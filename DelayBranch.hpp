 /* 
 * * Description:
 * * This file contains a delay line class used to represent a sample
 * * of sound travelling through a waveguide in either fwd (S+) or bwd (S-)
 * * direction. Two DelayBranch objects in opposite directions form one
 * * bidirectional tube/string.
 * *
 * * Author:
 * * JEP J. Enrique Peraza, P&G Labs, LLC
 * *
 */
#pragma once
#include "DelayLine.hpp"
#include "FarrowDelayLine.hpp"
#include "KernelDeinterpolator.hpp"
#include "WGTypes.hpp"

namespace sig::wg
{
  template<size_t MaxLen = 1 << 15>
  class DelayBranch final:public Node
  {
    public:
    DelayBranch(void) noexcept
    {
      this->dl=new sig::DelayLine<wg::Sample<float>,MaxLen>(); // Create a new delay line object.
      farrow=new FarrowDelayLine<float,MaxLen,3>(); // Create a new Farrow delay line object.
      assert(this->dl && "DelayBranch: Failed to allocate delay line!"); // Ensure the delay line was allocated successfully.
      this->mu=0.0f; // Initialize the fractional delay to 0.
    }
    ~DelayBranch(void) noexcept
    {
      delete this->dl; // Delete the delay line object.
      this->dl=nullptr; // Set the delay line pointer to null.
      delete farrow; // Delete the Farrow delay line object.
      farrow=nullptr; // Set the Farrow delay line pointer to null.
    }
    bool SetMu(float m) noexcept
    {                                   // ----------- SetMu ------------- //
      if (m<0.5f||m>1.5f) return false; // Ensure the fractional delay is in range.
      this->mu=m;                       // Set the fractional delay.
      return true;                      // Return true if successful.
    }                                   // ----------- SetMu ------------- //
    bool SetDelay(size_t d) noexcept { return (this->dl->SetDelay(d)); }
    float Read(void) const noexcept
    {
      return farrow->ReadFrac(mu); // Read a sample from the delay line with fractional delay.
    }
    // Write a sample to tail of the delay line:
    inline void Write(wg::Sample<float> s) noexcept { WriteFrac(0,s); }
    inline void WriteFrac(
      size_t widx=0,
      double s=0.0,
      float mu=0.0f) noexcept
    {                                   // ----------- WriteFrac ------------- //
      if (dl!=nullptr)
        KernelDeinterpolator<float,MaxLen>::Deposit(
          dl,widx,s,mu); // Write the sample to the delay line with interpolation.
    }                                   // ----------- WriteFrac ------------- //
    // Propagate the delay line by 'n' samples.
    // repeatedly read the head and write it back to the tail.
    void Propagate(size_t n) noexcept override
    {                                   // ----------- Propagate ----------- //
        for (size_t i=0;i<n;++i)        // For each sample to propagate...
        {
          //this->dl->Write(this->dl->Read());// Propagate the sample across the line.
          float s=this->dl->Read(); // Read the sample from the delay line.
          this->dl->Write(s);       // Write the sample back to the delay line.
        }  
    }                                   // ----------- Propagate ----------- //
    
    inline size_t GetDelay(void) const noexcept { return this->dl->GetDelay(); }
    inline size_t GetMaxLen(void) const noexcept { return MaxLen; }
    float Peek(void) const noexcept
    {
      return farrow->Read();
    }
    void Clear(void) noexcept
    {
      this->dl->Clear(); // Clear the delay line.
    }                                   // Clear the delay line.
    private:
      DelayLine<float,MaxLen>* dl;
      FarrowDelayLine<float,MaxLen,3>* farrow;// The delay line object.
      float mu{0.3f}; // The fractional delay in samples.
      //static_assert(MaxLen > 0 && (MaxLen & (MaxLen - 1)) == 0, "MaxLen must be a power of two for efficient wrap-around.");
      //double len{MaxLen};             // The length of the delay line in samples.
   };
} // namespace sig::wg
