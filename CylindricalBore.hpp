 /* 
 * * P&G Labs, LLC
 * *
 * * Filename: CylindricalBore.hpp
 * * Description:
 * * This file contains a CylindricalBore class that represents a cylindrical bore
 * * in a waveguide network.
 * * It is used to model the acoustic properties of a cylindrical bore in a waveguide.
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <cassert>
#include <algorithm>
#include <cmath>
#include "DelayBranch.hpp"
#include "FilterFactory.hpp"
namespace sig::wg {
template <size_t MaxLen = 1 << 15>
  class CylindricalElement final: public Node
  {
    public: 
      CylindricalElement(void) noexcept
      {
        fwd=new DelayBranch<MaxLen>;        // Create a new forward delay branch.
        rev=new DelayBranch<MaxLen>;        // Create a new reverse delay branch.
        damp=FilterFactory<Sample<float>>::ButterWorthLP(fs,fc,Q); // Create a damping filter with the specified parameters.
        fs=48000.0; // Default sampling frequency.
        f0=440.0; // Default fundamental frequency.
        fc=22050.0; // Default cutoff frequency.
        length=1.0; // Default length of the bore in meters.
        radius=0.1; // Default radius of the bore in meters.
        press=0; // Initialize the pressure to 0.
      }
      ~CylindricalElement(void) noexcept override
      {
        delete fwd; // Delete the forward delay branch.
        delete rev; // Delete the reverse delay branch.
        fwd=nullptr; // Set the forward delay branch pointer to null.
        rev=nullptr; // Set the reverse delay branch pointer to null.
      }                                   // Destructor for CylindricalElement.
    bool Prepare(double fs, double length, double radius)
    {                                   // ----------- Prepare ----------------- //
      
      this->fs=fs; // Set the sampling frequency.
      this->length=length; // Set the length of the bore.
      this->radius=radius; // Set the radius of the bore.
      // Number of samples for the length of the bore.
      size_t N=static_cast<size_t>(fs*length/c);
      if (N<1) N=1;                     // Ensure minimum length of 1 sample.
      fwd->SetDelay(N);                  // Set the delay for the forward branch.
      rev->SetDelay(N);                  // Set the delay for the reverse branch.
      // Pick a cutoff fractionally below the Nyquist frequency.
      double nyquist=fs*0.5;            // Nyquist frequency.
      this->fc=std::min(fc,nyquist*0.99);
      damp=FilterFactory<Sample<float>>::ButterWorthLP(fs,fc,this->Q);
      press=0;                     // Initialize output to 0.
      return true;                      // Return true if prepared successfully.
    }                                   // ----------- Prepare ----------------- //
    // Propagate the samples through the cylindrical bore.
    void Propagate(size_t /* n (unused)*/) noexcept override
    {                                   // ----------- Propagate ----------------- //
      auto infwd=rev->Read();           // Left going meets
      auto inrev=fwd->Read();           // Right going meets
      auto outF=this->damp.ProcessSample(-infwd);
      auto outR=this->damp.ProcessSample(-inrev);
      this->fwd->Write(outF);            // Write the processed sample back to the forward branch.
      this->rev->Write(outR);            // Write the processed sample back to the reverse branch.
      this->press=outF+outR;            // Sum of waves = physical displacement.
    }                                   // ----------- Propagate ----------------- //
    void Excite(Sample<float> s) noexcept { fwd->Write(s); rev->Write(-s); } // Write excitation to both branches.
    // Output the current position sample (output).
    Sample<float> Pressure(void) const noexcept { return press; } // Get the current position sample.
    // Getters for the forward and reverse branches.
    inline const DelayBranch<MaxLen>* GetFwd(void) const noexcept { return fwd; }
    inline const DelayBranch<MaxLen>* GetRev(void) noexcept { return rev; }
    inline const DelayBranch<MaxLen>*  const GetFwd(void) noexcept { return fwd; }
    inline const DelayBranch<MaxLen>* const GetRev(void) const noexcept { return rev; } // Get the reverse branch.
    // Get the damping filter.
    inline const BiQuad<Sample<float>>& GetDamp(void) const noexcept { return damp; } // Get the damping filter.
    inline void SetDamp(const BiQuad<Sample<float>>& d) noexcept { damp=d; } // Set the damping filter.
    inline const double GetSamplingFrequency(void) const noexcept { return fs; } // Get the sampling frequency.
    inline void SetSamplingFrequency(double fs) noexcept { this->fs=fs; } // Set the sampling frequency.
    inline const double GetFundamentalFrequency(void) const noexcept { return f0; } // Get the fundamental frequency.
    inline void SetFundamentalFrequency(double f0) noexcept { this->f0=f0; } // Set the fundamental frequency.
    inline const double GetCutoffFrequency(void) const noexcept { return fc; } // Get the cutoff frequency.
    inline void SetCutoffFrequency(double fc) noexcept 
    {
       this->fc=fc;
       damp=FilterFactory<Sample<float>>::ButterWorthLP(fs,fc,Q); // Update the damping filter with the new cutoff frequency. 
    }
    // Set the length of the bore.
    inline void SetLength(double length) noexcept { this->length=length; } // Set the length of the bore.
    // Get the length of the bore.
    inline double GetLength(void) const noexcept { return length; } // Get the length of the bore.
    // Set the radius of the bore.
    inline void SetRadius(double radius) noexcept { this->radius=radius; } // Set the radius of the bore.
    // Get the radius of the bore.
    inline double GetRadius(void) const noexcept { return radius; } // Get the radius of the bore.
    private:
      DelayBranch<MaxLen>* fwd; // Forward delay branch.
      DelayBranch<MaxLen>* rev; // Reverse delay branch.
      BiQuad<Sample<float>> damp; // Damping filter.
      Sample<float> press;        // Current pressure sample.
      const double c{343.0}; // Speed of sound in air (m/s).
      double f0{440.0}; // Fundamental frequency.
      double fs{48000.0}; // Sampling frequency.
      double fc{22050.0}; // Cutoff frequency.
      double length{1.0}; // Length of the bore (meters).
      double radius{1.0}; // Radius of the bore (meters).
      double Q{0.70710678}; // Quality factor for damping.
  };
}  // namespace sig::wg
