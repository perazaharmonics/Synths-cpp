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
pragma once
#include <cstddef>
#include <array>
#include "WGTypes.hpp"
#include "DelayBranch.hpp"
#include "FilterFactory.hpp"
#include "BiQuad.hpp"
namespace sig::wg
{
  template<typename T=float,size_t N,size_t D, T mut, T muf,size_t order=5>
  struct ScatteringJunction final:public Node
  {
    static_assert(N>=2&&N<=8,"juction size unreasonable.");
    // Attach N travelling-wave branch pointers at construction
    explicit ScatteringJunction(
      const std::array<DelayBranch<>*, N>& branches) noexcept
    : bra(branches) {}
    
    // No-op Prepare. Branch length and delays should be set externally.
    bool Prepare(double fs,double /*f0*/,double fc,size_t idelay, double m0, double m1)) noexcept
    {
      this->fs=fs; // Set the sampling frequency.
      this->fc=fc; // Set the cutoff frequency for the damping filter.
      mut=m0;       // Set Thrian Fractional delay
      muf=m1;       // Set Farrow fractional delay
      damp=ff.Bessel(fs,fc,0.7071f); // Create the damping filter.
      damp.Clear();  // Clear filter state so it doesn’t attenuate the signal.
      for (size_t i=0;i<N;++i) // For al branches
      {                        // configure fractional delays.
        bra[i].Prepare(idelay,mut,muf);// Prepare the waveguide
        bra[i].SetFractionalDelay(mut);// Thiran's fractional delay
        bra[i].SetMuFarrow(muf);// Farrow's Fractional Delay
      }                         // Done setting up the Fractional delays.
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
    size_t o{order};  // Order of interplation filters
    double fs{48000.0}; // Sampling frequency.
    double fc{22050.0}; // Cutoff frequency for the damping filter.
  }; // end of ScatteringJunction
} // namespace sig::wg