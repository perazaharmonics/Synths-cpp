 /* 
 * * P&G Labs, LLC
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
  template<size_t N>
  struct ScatteringJunction final:public Node
  {
    static_assert(N>=2&&N<=8,"juction size unreasonable.");
    // Attach N travelling-wave branch pointers at construction
    explicit ScatteringJunction(
      const std::array<DelayBranch<>*, N>& branches) noexcept
    : bra(branches) {}
    
    // No-op Prepare. Branch length and delays should be set externally.
    bool Prepare(double fs,double /*f0*/,double fc) noexcept
    {
      this->fs=fs; // Set the sampling frequency.
      this->fc=fc; // Set the cutoff frequency for the damping filter.
      damp=ff.Bessel(fs,fc,0.7071f); // Create the damping filter.
      damp.Clear();  // Clear filter state so it doesn’t attenuate the signal.
      return true; // Return true to indicate successful preparation.
    } 
    
    // -------------------------------- //
    // Ideal junction reflection rule: S+[k]<---2*rho-S-[k], where
    // rho is the average of incident pressure waves (S-) The scattering matrix.
    // -------------------------------- //
    void Propagate(size_t n) noexcept override
    {                                   // -------- Propagate --------------- //
        //for (size_t i=0;i<N;i++)
        //  bra[i]->Propagate(n);          // Propagate each branch by 'n' samples.  
      // Propagate each branch by 'n' samples.
        // ---------------------------- //
        // Read the current head samples (incident waves)
        // ---------------------------- //
        std::array<float,N> inc{};      // Incident pressure waves (S-).
        float sum=0.0f;                 // Sum of incident pressure waves.
        for (size_t i=0;i<N;++i)         // For each branch...
        {
          inc[i]=bra[i]->Peek();          // Read the incident pressure wave.
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
          bra[i]->Write(out);             // Write the scattered pressure wave to the branch.
        }                                // Done writing scattered pressure waves.
        // Advance each branch by n samples again.
        for (size_t i=0;i<N;++i)
          bra[i]->Propagate(n);          // Propagate each branch by 'n' samples again.
    }                                   // -------- Propagate --------------- //
  private:
    std::array<DelayBranch<>*,N> bra;      // The branches connected to the junction.
    FilterFactory<float> ff; // Filter factory for creating filters.
    BiQuad<float> damp; // Damping filter for the junction.
    double fs{48000.0}; // Sampling frequency.
    double fc{22050.0}; // Cutoff frequency for the damping filter.
  }; // end of ScatteringJunction
} // namespace sig::wg
