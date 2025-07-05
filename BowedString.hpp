 /* 
 * * 
 * *
 * * Filename: BowedString.hpp
 * * Description:
 * * This file contains a BowedString class that represents a bowed string
 * * in a waveguide network.
 * * It is used to model the acoustic properties of a bowed string in a waveguide.
 * * It uses the Smith-Hausser Elasto-plastic friction model where:
 * * sigma = bow position; vrel = v_string - v_bow;
 * * v_bow = bow velocity; v_string = string velocity.
 * * The friction force is computed as:
 * * F= Fb & tanh(vrel/v0) (Good enough for layer-2)
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <cassert>
#include <cmath>
#include "WGTypes.hpp"
#include "FarrowDelayLine.hpp"
#include "OnePole.hpp"
namespace sig::wg
{
  struct BowedString :public Node
  {
    static constexpr size_t MaxLen=2048; // Must be power of two
    // Reset the state of the bowed string.
    // Prepare the bowed string with sampling frequency, fundamental frequency, and loop cutoff frequency.
    bool Prepare(double fs, double f0,double pos, double loopfc,double loopQ) noexcept
    {                                  // ----------- Prepare ----------------- //
        this->fs=fs;                   // Set the sampling frequency.
        this->f0=f0;                   // Set the fundamental frequency.
        this->loopFc=loopfc;           // Set the loop cutoff frequency.
        this->loopQ=loopQ;             // Set the loop quality factor.
        this->pos=pos;                 // Set the bow position (0.0 to 1.0).
        // Clear the reverse delay line.
        size_t N=static_cast<size_t>(fs/(2.0*f0)); // Half period of the string element.
        // Assert half period
        assert(N < MaxLen && "MaxLen too small for this f0");
        double Ni=fs/(2.0*f0); // Half period of the string element in samples.
        this->revdel.SetDelay(static_cast<size_t>(std::round(Ni)));      // Set the reverse delay line length.
        this->loss.SetConf(OnePole<float>::Conf::Lowpass);// Set the loss filter to lowpass.
        // --------------------------- //
        // Bridge <--> Bow <---> nut
        // --------------------------- //
        this->fwddel.RampTo(pos*N,1);  // Ramp the forward delay line to the bow position.
        this->revdel.RampTo((1.0f-pos)*N,1); // Ramp the reverse delay line to the nut position.
        this->loss.Prepare(fs,loopfc); // Prepare the loss filter with the sampling frequency and loop cutoff frequency.
        return true;                   // Return true if prepared successfully.
    }                                  // ----------- Prepare ----------------- //
    void ResetState(void) noexcept
    {
      this->fwddel.Clear();            // Clear the forward delay line.
      this->revdel.Clear();            // Clear the reverse delay line.
      this->loss.Reset();              // Reset the loss filter state.
      this->vBow=0.1f;                 // Reset the bow velocity.
      this->Pm=0.3f;                   // Reset the mouth pressure.
      this->Pmax=10.f;                 // Reset the maximum pressure.
      this->S=1e-4f;                   // Reset the bow area.
      this->mu=1.f;                    // Reset the slope parameter for the friction model.
      this->rho=1.21f;                 // Reset the density of the string.
      this->c=343.f;                   // Reset the speed of sound in the string.
    }                                  // ----------- ResetState ----------------- //
    // Per-sample nonlinear scattering.
    void Propagate(size_t n) noexcept override
    {
      for (size_t i=0;i<n;++i)          // For each sample to propagate...
      {
        //  compute string velocity at bow (before writing new samples.
        float vstring=(fwddel.Read()+revdel.Read());// Velocity proxy
        float vRel=vstring-vBow;                     // Relative velocity between string and bow.
        float dP=std::max(0.f, std::min(float(Pm-vstring), float(Pmax)));
        float U=(dP>0&&dP<Pmax)?S*mu*dP*std::pow(1-dP/Pmax,2):0.f; // Volume flow based on pressure difference.
        if (Pmax <= 0.f) Pmax = 1.f;                // avoid div-by-zero
        if (Pm >= Pmax) Pmax = Pm * 1.1f;           // guarantee dP < Pmax

        // ----------------------------------------- //
        // Convert to pressure.
        // ----------------------------------------- //
        float Pw=U*(rho*c/S);                        // Convert volume flow to pressure wave using characteristic impedance.
        // Inject friction in both travelling waves.
        float inR=revdel.Read()+Pw;                  // Read the reverse delay line and add the pressure wave.
        float inF=fwddel.Read()-Pw;                  // Read the forward delay line and subtract the pressure wave.
        fwddel.Write(loss.ProcessSample(inR));       // Write the processed sample back to the forward delay line.
        revdel.Write(loss.ProcessSample(inF));       // Write the processed sample back to the reverse delay line.
        // ----------------------------------------- //
        // Update output and fractional read delay pointers.
        // ----------------------------------------- //
        this->fwddel.Tick();                   // Update the forward delay line.
        this->revdel.Tick();                   // Update the reverse delay line.
        // ----------------------------------- //
        // Update the output sample.
        // ----------------------------------- //
        this->out=fwddel.Read()+revdel.Read(); // Update the output sample as the sum of both delay lines.
      }                                // Done processing the samples.
    }                                   // ----------- Propagate ----------------- //
    float Pressure(void) const noexcept { return this->out; } // Get the current pressure sample (output).
    void SetMouthPressure(float v) noexcept { this->vBow=v; } // Set the bow velocity.
    void SetBowVelocity(float v) noexcept { this->vBow=v; } // Set the bow velocity.
    void SetBowPosition(float p) noexcept
    {
      this->pos=p;                   // Set the bow position (0.0 to 1.0).
      size_t N=static_cast<size_t>(fs/(2.0*f0)); // Half period of the string element.
      fwddel.RampTo(p*N,1); // Ramp the forward delay line to the new bow position.
      revdel.RampTo((1.0f-p)*N,1); // Ramp the reverse delay line to the new nut position.
    }
    // Loop alpha API.
    void SetLoopCutoff(double fc) noexcept {loopFc = fc; loss.Prepare(fs, fc); } // Set the loop cutoff frequency.
    double GetLoopCutoff(void) const noexcept { return this->loopFc; } // Get the loop cutoff frequency.
    void SetLoopQ(double q) noexcept {  } // Set the loop quality factor.
    double GetLoopQ(void) const noexcept { return this->loopQ; } // Get the loop quality factor.
    // Friction API.
    void SetMaxFriction(float m) noexcept { Pmax=m; } // Set the maximum friction force (scaled).
    void SetBernoulli(float m) noexcept { mu=m; } // Set the slope parameter for the friction model.
    void SetArea(float a) noexcept { this->S=a; } // Set the bow area (not used in this model).
    void SetDensity(float d) noexcept { this->rho=d; } // Set the density of the string (not used in this model).
    void SetSpeed(float s) noexcept { this->c=s; } // Set the speed of sound in the string (not used in this model).
    private:
      double fs=44100.0;                 // Sampling frequency (Hz).
      double f0=440.0;                   // Fundamental frequency (Hz).
      double loopFc=5000.0;             // Loop cutoff frequency (Hz).
      double loopQ=0.707;                // Loop quality factor.
      FarrowDelayLine<float,1024UL> fwddel,revdel; // Forward and reverse delay lines.
      OnePole<float> loss;               // Lossy filter for the string element.
      float vBow=0.1f;                // Bow velocity (m/s).
      float Pm=0.3f;                  // Mouth pressure (Pa).
      float Pmax=10.f;              // Maximum pressure (Pa).
      float S=1e-4f;               // Bow area (m^2).
      float mu=1.f;                // Slope parameter for the friction model.
      float rho=1.21f;             // Density of the string (kg/m^3).
      float pos=0.1f;               // Bow position (0.0 to 1.0).
      float c=343.f;               // Speed of sound in the string (m/s).
      float out=0.f;               // Output pressure sample.
    };
}
