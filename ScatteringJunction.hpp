/* 
 * *
 * * Filename: ScatteringJunction.hpp
 * * Description:
 * * This file contains a scattering juction class. The scattering junction
 * * is an abrtraction of a waveguide junction, where multiple waveguides
 * * can connect and scatter energy through. They are either all-pass
 * * or lossy, depending on the configuration. These are 2D  3-port T junctions
 * * or 4-port cross junsctions. Although its templated for N ports,
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <cstddef>
#include <array>
#include "WGTypes.hpp"
#include "DelayBranch.hpp"
#include "FilterFactory.hpp"
#include "BiQuad.hpp"
namespace sig::wg
{
  template<
  typename T=double,                     // Data type we are processing
  size_t N,                             // Number of branches for the scattering junction
  size_t MaxLen=1<<15,                  // Maximum length of the delay line in samples
  size_t K=5,                           // Thiran filter order
  size_t P=5>                           // Farrow filter order
  struct ScatteringJunction final:public Node
  {
    static_assert(N>=2&&N<=8,"juction size unreasonable.");
    // Attach N travelling-wave branch pointers at construction
    explicit ScatteringJunction(
      const std::array<DelayBranch<>*, N>& branches) noexcept
    : bra(branches) {}
    
    // Set the Scattering junction / waveguide parameters,and user settings. Then
    // prime the delay branches with zeroes. 
    bool Prepare(double fs,double /*f0*/,double fc,size_t idelay, double m0, double m1) noexcept
    {
     if (fs<=0.0||fc>=0.5*fs||P!=K||idelay<K) return false; // Sanitize input!
      this->fs=fs; // Set the sampling frequency.
      this->fc=fc; // Set the cutoff frequency for the damping filter.
      idelay=std::max<size_t>(1,idelay); // Sanitize the integer delay.
      mut=m0;       // Set Thrian Fractional delay
      muf=m1;       // Set Farrow fractional delay
      damp=ff.Bessel(fs,fc,0.7071f); // Create the damping filter.
      damp.Clear();  // Clear filter state so it doesn?t attenuate the signal.
      for (size_t i=0;i<N;++i) // For al branches
        bra[i]->Prepare(idelay,mut,muf);// Prepare the waveguide
      // Prime waveguides according to Group Delay
      size_t maxlat=idelay+static_cast<size_t>(muf+mut);
      for (size_t i=0;i<maxlat;++i) // For the max group delay
        Propagate(T(0)); // Circulate an impulse
      return true; // Return true to indicate successful preparation.
    } 
    
    // -------------------------------- //
    // Ideal junction reflection rule: S+[k]<---2*rho-S-[k], where
    // rho is the average of incident pressure waves (S-) The scattering matrix.
    // -------------------------------- //
    void Propagate(size_t n) noexcept override
    {                                   // -------- Propagate --------------- //
        // ---------------------------- //
        // Read the current head samples (incident waves)
        // ---------------------------- //
        std::array<float,N> inc{};      // Incident pressure waves (S-).
        float sum=0.0f;                 // Sum of incident pressure waves.
        for (size_t i=0;i<N;++i)         // For each branch...
        {
          inc[i]=bra[i]->Read();          // Use Read instead of Peek to read the incident pressure wave.
          sum+=inc[i];                   // Sum the incident pressure waves.
        }                                // Done reading incident pressure waves.
        // ---------------------------- //
        // Calculate the average incident pressure wave (rho).
        // ---------------------------- //
        float rho=sum/static_cast<float>(N); // Average incident pressure wave.
        // ---------------------------- //
        // Scatter into the tail of each line
        // ---------------------------- //
        for (size_t i=0;i<N;++i)        // For each branch...
        {
          float out=rho*2.f-inc[i];
          out=damp.ProcessSample(out); // Apply damping to the output.
          bra[i]->Write(out);             // Use Write instead of WriteFrac to write the scattered pressure wave to the branch.
        }                                // Done writing scattered pressure waves.
        // Advance each branch by n samples again.
        for (size_t i=0;i<N;++i)
          bra[i]->Propagate(n);          // Propagate each branch by 'n' samples again.
    }                                   // -------- Propagate --------------- //
  private:
    std::array<DelayBranch<>*,N> bra;      // The branches connected to the junction.
    FilterFactory<float> ff; // Filter factory for creating filters.
    BiQuad<float> damp; // Damping filter for the junction.
    double mut{0.25}; // Fractional delay interpolation factor for Thiran
    double muf{0.48}; // Fractional delay for Farrow
    size_t idelay{5}; // The integer delay
    size_t o{K};  // Order of interplation filters
    double fs{48000.0}; // Sampling frequency.
    double fc{22050.0}; // Cutoff frequency for the damping filter.
  }; // end of ScatteringJunction
} // namespace sig::wg
