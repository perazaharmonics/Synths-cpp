/*
* *
* * Filename: Oscillator.hpp
* * 
* * Description: Musical oscillators for sound synthesis. It contains the following base oscillators:
* *  - SineOscillator: A simple sine wave oscillator.
* *  - SawOscillator: A bandlimited sawtooth oscillator using PolyBLEP
* *  - SquareOscillator: A square wave oscillator using PolyBLEP
* *  - TriangleOscillator: A triangle wave oscillator is just the integral of a square wave.
* *  - WaveTableOscillator: A wave table oscillator that can generate any waveform from a precomputed table.
* *  - MakeWaveTable: A utility function to create a sine wave table for a given table size.
* *
* * Author: JEP J. Enrique Peraza
* *
*/
#pragma once
#include <cmath>
#include <cstdint>
#include <vector>
namespace sig::osc
{
    // -------------------------------- //
    // Musical oscillator base class
    // -------------------------------- //
    template <typename T=float>
    class Oscillator
    {
      public:
        Oscillator(void):
            fs(48000.0), // Default sampling frequency
            f0(440.0),   // Default fundamental frequency
            phase(0.0),  // Initial phase
            phinc(0.0)   // Phase increment
        {                                   // Default constructor
            UpdatePhaseInc();               // Update the phase increment
        }                                   // Default constructor
        virtual ~Oscillator(void) noexcept = default; // Default destructor
        // Set the sampling frequency and fundamental frequency
        void SetSamplingFrequency(double fs) noexcept { this->fs=fs; UpdatePhaseInc(); }
        void SetFundamentalFrequency(double f0) noexcept { this->f0=f0; UpdatePhaseInc(); }
        void SetAmplitude(double a) noexcept { this->amp=a; }
        void ResetPhase(void) noexcept { this->phase=0.0; } // Reset the phase to zero
        /// Return one sample of the oscillator output
        virtual T Process(void) noexcept
        {
            phase+=phinc;          // Increment the phase
            if (phase>=twoPi)      // If phase exceeds 2p
            phase-=twoPi;        // Wrap around the phase
            return static_cast<T>(0);
        }
      protected:
        double fs{48000.0}; // Sampling frequency
        double f0{440.0}; // Fundamental frequency
        double phase{0.0}; // Current phase
        double phase{0.0}; // Current phase
        double phinc{0.0}; // Phase increment per sample
        double amp{1.0};  // Amplitude of the oscillator
        static constexpr double twoPi=double(2*M_PI);
        virtual void UpdatePhaseInc(void) noexcept
        {                                   // Update the phase increment based on the sampling frequency and fundamental frequency
            phinc=twoPi*f0 /fs;     // Calculate the phase increment
        }                                   // Update the phase increment
    };
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  // Sine Oscillator
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  template <typename T=float>
  class SineOscillator: public Oscillator<T>
  {
    public:
      using Osc=Oscillator<T>::Oscillator; // Inherit from Oscillator
      T Process(void) noexcept override
      {                                 // ~~~~~~~~~ Process ~~~~~~~~~~ //
        this->phase+=this->phaseinc;    // Increment the phase
        if (this->phase>=Osc::twoPi)    // If phase exceeds
          this->phase-=Osc::twoPi;      // Wrap around the phase
        return this->amp*std::sin(this->phase); // Return sine wave.
      }                                 // ~~~~~~~~~ Process ~~~~~~~~~~ //
  }; // Sine Oscillator
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  // Sine Wave Table Generator
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  template <typename T=float>
  inline std::vector<T> MakeWaveTable (
    size_t s=2048) noexcept               // The size of the sine wave table 
  {                                       // ~~~~~~ GenerateSineWaveTable ~~~~~~ //
    std::vector<T> w(s);                  // Create a vector of size s
    for (size_t i=0;i<s;++i)              // For each entry in the table
      w[i]=std::sin(twoPi*T(i)/T(s));     // Calculate the sine value
    return w;                             // Return the completed table
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  // Bandlimited Saw using PolyBLEP
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  template <typename T=float>
  class SawOscillator:public Oscillator<T>
  {
    public:
      using Osc=Oscillator<T>::Oscillator;
      SawOscillator(void) noexcept { }
      T Process(void) noexcept override
      {                                 // ~~~~~~~~~ Process ~~~~~~~~~~ //
        phase+=phinc;                   // Increment the phase
        if (phase>=Osc::twoPi)          // If phase exceeds 2pi
          phase-=Osc::twoPi;            // Wrap around the phase
        T t=phase/Osc::twoPi;           // Normalize the phase to [0,1]
        T dt=phinc/Osc::twoPi;          // Calculate the delta time
        T y=T(2)*t-T(1);                // Sawtooth wave: y = 2t - 1
        if (t<dt)                       // If t is less than dt
          y-=PolyBlep<T>(t,dt);         // Apply PolyBLEP for
        else if (t>T(1)-dt)             // If t is greater than 1-dt
          y-=PolyBlep<T>(t-T(1),dt);    // Apply PolyBLEP for the falling edge
        return amp*y;                   // Return the processed sample
      }                                 // ~~~~~~~~~ Process ~~~~~~~~~~ //
    private:
      template <typename U>
      static U PolyBlep(                // PolyBLEP function for bandlimited sawtooth wave
        U t,                            // time in seconds
        U dt)                           // Time step or delta time
        {                               // ~~~~~~~~~ PolyBlep ~~~~~~~~~~ //
          if (t<dt)                     // Is the time less than delta time?
          {                             // Yes
            U x=t/dt;                   // Calculate the normalized time
            return x+x=x*x-U(1);        // Return the PolyBLEP value
          }                             // 
          else if (t>U(1)-dt)           // Is the time greater than 1-delta time?
          {                             // Yes
            U x=(t-U(1))/dt;            // Calculate the normalized time
            return x*x+x+x+U(1);        // Return the PolyBLEP value
          }                             // Otherwise
          return U(0);                  // Return zero
        }                               // ~~~~~~~~~ PolyBlep ~~~~~~~~~~ //
  }; // Saw Oscillator
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  // Square Oscillator
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  template <typename T=float>
  class SquareOscillator: public Oscillator<T>
  {
    public:
      using Osc=Oscillator<T>::Oscillator;
      SquareOscillator(void) noexcept { }
      void SetDutyCycle(T d) noexcept { duty=std::clamp(d,T(0.001),T(0.999)); } // Set the duty cycle of the square wave
      T Process (void) noexcept override
      {                                 // ~~~~~~~~~~~~~ Process ~~~~~~~~~~~~~~ //
        phase+=phinc;                   // Increment the phase
        if (phase>=Osc::twoPi)          // If phase exceeds 2pi
          phase-=Osc::twoPi;            // Wrap around the phase
        T t=phase/Osc::twoPi;           // Normalize the phase to [0,1]
        T dt=phinc/Osc::twoPi;          // Calculate the delta
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Get the rising and falling edges of the square wave
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        T y=(t<duty)?T(1):T(-1);        // Square wave: 1 for t < duty, -1 for t >= duty
        if (t<dt)                       // If t is less than dt
          y+=PolyBlep<T>(t,dt);         // Apply PolyBLEP for the rising edge
        else if (t<duty&&t>duty-dt)     // Is the in the falling edge region?
          y-=PolyBlep<T>(t-duty,dt);    // Apply PolyBLEP for the falling edge
        else if (t>=duty&&t<duty+dt)    // Is the in the rising edge region?
          y+=PolyBlep<T>(t-duty,dt);    // Apply PolyBLEP for the rising edge
        return amp*y;                   // Return the processed sample
      }                                 // ~~~~~~~~~~~~~ Process ~~~~~~~~~~~~~~ //
    private:
      T duty{0.5};                      // The duty cycle of the square wave, default is 0.5
      template <typename U>
      static U PolyBlep (               // PolyBLEP function for bandlimited square wave
        U t,                            // time in seconds
        U dt) noexcept                  // Time step or delta time    
      {                                 // ~~~~~~~~~ PolyBlep ~~~~~~~~~~ //
        if (t<dt)                       // Is the time less than delta time?
        {                               // Yes
          U x=t/dt;                     // Calculate the normalized time
          return x*x+x+x-U(1);          // Return the PolyBLEP value
        }                               // 
        else if (t>U(1)-dt)             // Is the time greater than 1-delta time?
        {                               // Yes
          U x=(t-U(1))/dt;              // Calculate the normalized time
          return x*x+x+x+U(1);          // Return the PolyBLEP value
        }                               // Otherwise
        return U(0);                    // Return zero
      }                                 // ~~~~~~~~~ PolyBlep ~~~~~~~~~~ //
  }; // Square Oscillator
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  // Triangle Oscillator
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  template <typename T=float>
  class TriangleOscillator: public Oscillator<T>
  {
    using Osc=Oscillator<T>::Oscillator; // Inherit from Oscillator
    public:
      TriangleOscillator(void) noexcept { }
      T Process(void) noexcept override
      {                                 // ~~~~~~~~~~ Process ~~~~~~~~~~ //
        prev+=phinc;                    // Store phase increment in previous
        if (prev>=Osc::twoPi)           // If phase exceeds 2pi
          prev-=Osc::twoPi;             // Wrap around the phase
        T sq=PolySquare<T>(prev);       // Get the square wave value.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // The integral of the square wave is a triangle wave
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        inte+=sq*(phinc/Osc::twoPi);    // Integrate the square wave
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Wrap integrator to [-1,1]
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        if (inte>T(1)) inte-=T(2);      // If integrator is greater than 1
        else if (inte<T(-1)) inte+=T(2);// If integrator is less than -1
        return amp*inte;                // Return the triangle wave
      }                                 // ~~~~~~~~~~ Process ~~~~~~~~~~ //
    private:
      T prev=T(0);                      // Previous phase increment.
      T inte=T(0);                      // Integrator
      T PolySquare(T ph) noexcept       // Takes the phase
      {                                 // ~~~~~~~~~ PolySquare ~~~~~~~~~~ //
        // Same logic as square without amplitude....
        T y=ph/Osc::twoPi<T(0.5)?T(1):T(-1);
        return y;                       // Return the square wave value with no edge corrections
      }                                 // ~~~~~~~~~ PolySquare ~~~~~~~~~~ //
  }; // Triangle Oscillator
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  // Wavetable Oscilator
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
  template <typename T=float>
  class WavetableOscillator: public Oscillator<T>
  {
    using Osc=Oscillator<T>::Oscillator; // Inherit from Oscillator
    public:
      WavetableOscillator(const std::vector<T>& tb={}):tbl(tb) { }
      void SetWavetable(const std::vector<T>& tb) noexcept { tbl=tb; } // Set the wavetable
      T Process(void) noexcept override
      {                                 // ~~~~~~~~~ Process ~~~~~~~~~~ //
        if (tbl.empty()) return T(0);   // If the wavetable is empty, return zero
        phase+=phinc;                   // Increment the phase
        if (phase>=Osc::twoPi)          // If phase exceeds 2pi
          phase-=Osc::twoPi;            // Wrap around the phase
        T idx=phase/Osc::twoPi*T(tbl.size());// Index into the wavetable
        size_t i0=static_cast<size_t>(std::floor(idx))%tbl.size(); // Get the integer part of the index
        size_t i1=(i0+1)%tbl.size();    // Get the next index in the wavetable
        T f=idx-std::floor(idx);        // Get the fractional part of the index
        return amp*(tbl[i0]+(tbl[i1]-tbl[i0])*f);// Linear interpolation between the two samples
      }                                 // ~~~~~~~~~ Process ~~~~~~~~~~ //
    private:
      std::vector<T> tbl;               // The wavetable
  }; // Wavetable Oscillator
} // namespace sig::osc
