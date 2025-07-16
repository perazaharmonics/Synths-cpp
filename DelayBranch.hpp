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
#include <algorithm>
#include "DelayLine.hpp"
#include "FarrowDelayLine.hpp"
#include "WGTypes.hpp"

namespace sig::wg
{
  template<size_t MaxLen = 1 << 15>
  class DelayBranch final:public Node
  {
    public:
    DelayBranch(void) noexcept
    {
      farrow=new FarrowDelayLine<float,MaxLen,3>(); // Create a new Farrow delay line object.
      farrow->Prepare(mu);
    }
    ~DelayBranch(void) noexcept
    {
      delete farrow; // Delete the Farrow delay line object.
      farrow=nullptr; // Set the Farrow delay line pointer to null.
    }
    void Prepare(
      float d) noexcept
    {                                   // ----------- Prepare ------------- //
      farrow->Prepare(d); // Prepare the delay line for processing with the specified delay.
      this->mu=farrow->GetMu(); // Get the fractional delay from the Farrow delay line.
    }                                   // ----------- Prepare ------------- //
    bool SetMu(float m) noexcept
    {                                   // ----------- SetMu ------------- //
      farrow->SetMu(m); // Set the fractional delay in the Farrow delay line.
      return true;                      // Return true if successful.
    }                                   // ----------- SetMu ------------- //
    float GetMu(void) const noexcept
    {                                   // ----------- GetMu ------------- //
      return farrow->GetMu(); // Get the fractional delay from the Farrow delay line.
    }                                   // ----------- GetMu ------------- //
    float Read(void)
    {
      return farrow->Read();
    }
    float ReadFrac(float m) const noexcept
    {
      return farrow->ReadFrac(m); // Read a sample from the delay line with fractional delay.
    }
    // Write a sample to tail of the delay line:
    inline void Write(const wg::Sample<float>& s) noexcept 
    { 
      farrow->Write(s); // Write the sample to the delay line.
    }

    inline void WriteFrac(
      size_t widx=0,
      double s=0.0,
      float mu=0.0f) noexcept
    {                                   // ----------- WriteFrac ------------- //
      if (farrow!=nullptr)
        farrow->Write(static_cast<float>(s)); // Write the sample to the delay line.
    }                                   // ----------- WriteFrac ------------- //
    // Propagate the delay line by 'n' samples.
    // repeatedly read the head and write it back to the tail.
    void Propagate(size_t n) noexcept override
    {                                   // ----------- Propagate ----------- //
        farrow->Propagate(n); // Propagate the delay line by 'n' samples.
    }                                   // ----------- Propagate ----------- //
    float PeekIndex(size_t idx) const noexcept 
    {
      return farrow->PeekIndex(idx); // Peek at the sample at the specified index in the delay line.
    }
    inline size_t GetDelay(void) const noexcept { return farrow->GetDelay(); }
    inline void SetDelay(float f) noexcept
    {
      farrow->SetDelay(f);
    }
    float Peek(void) const noexcept
    {
      return farrow->Read();
    }
    void Clear(void) noexcept
    {
      farrow->Clear();
    }                                   // Clear the delay line.
    private:
      FarrowDelayLine<float,MaxLen,3>* farrow{nullptr}; // The delay line object.
      float mu{0.0f}; // The fractional delay in samples.
      //static_assert(MaxLen > 0 && (MaxLen & (MaxLen - 1)) == 0, "MaxLen must be a power of two for efficient wrap-around.");
      //double len{MaxLen};             // The length of the delay line in samples.
   };
} // namespace sig::wg
