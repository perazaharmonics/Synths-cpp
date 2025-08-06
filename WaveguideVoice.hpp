/*
* *
* * Filenaame: WaveguideVoice.hpp
* *
* * Description: WaveguideVoice is a waveguide-based voice synthesis engine
* * that simulates a string instrument using a waveguide model. It also serves
* * as an example of how to use an excitation signal to drive the waveguide.
* * The voice is excited by a pluck signal, and the waveguide simulates the
* * propagation of sound through a string.
* *
* * Author:
* * JEP J. Enrique Peraza
* *
*/
#pragma once
#include <cmath>
#include <memory>
#include <iostream>
#include "StringElement.hpp"
#include "BiQuad.hpp"
#include "FilterFactory.hpp"
#include "AudioBuffer.hpp"
#include "IProcessor.hpp"

namespace sig::wg
{
  enum class ToneType {Lowpass,Highpass,Bandpass,Shelfhigh,Shelflow};
  class WaveguideVoice
  {
    static constexpr size_t Maxlen=2048;       // Maximum length of the delay line in samples
    static constexpr size_t K=3;               // Farrow filter order
    static constexpr size_t P=3;               // Thiran filter order
    using String=StringElement<float,Maxlen,K,P>; // Our string element type
    public:
      // VoicePool Requirements
      inline int MidiNote(void) const noexcept { return cmidi; } // Get the current MIDI note id
      inline void NoteOff(void) noexcept { active=false; } // Turn off the voice
      inline bool IsActive(void) const noexcept { return active; } // Check if the voice is active
      bool Prepare(
        double fs,                          // Sample rate
        size_t b,                           // Block size
        size_t c) noexcept                  // Number of channels
      {                                     // ~~~~~~~~~ Prepare ~~~~~~~~~~ //
        if (fs<=0.0||b==0||c==0) return false; // Sanitize sample rate, block size and channels
        this->fs=fs;                        // Set the sample rate
        bs=b;                               // Set the block size
        ch=c;                               // Set the number of channels
        samples=bs*ch;                      // Set the number of samples to process
        string=std::make_unique<String>();  // Create a new StringElement instance
        if (!string->Prepare(fs,f0,loopfc,alpha,stage,D,mut,muf,order,pos))
          return false;                     // Preparation failed
        Assemble();                         // Assemble the tone shaping filter
        age=0;                              // Reset the age of the voice
        active=false;                       // Signify that the voice is not active
        return true;                        // Preparation successful
      }                                     // ~~~~~~~~~ Prepare ~~~~~~~~~~ //
      void ResetState(
        int mnote,                          // MIDI note id
        double mvel) noexcept               // The velocity of it.
      {
        cmidi=mnote;                        // Set midi note id.
        vel=mvel;                           // Set the velocity of the pluck
        amp=vel*scale;                      // Set the amplitude of the pluck
        f0=440.0*std::pow(2.0,(cmidi-69)/12.0); // Set the fundamental frequency
        string->Prepare(fs,f0,loopfc,alpha,stage,D,mut,muf,order,pos); // Prepare the string element
        string->Excite(amp);                // Excite the string with the pluck signal
        age=0;                              // Reset the age of the voice
        active=true;                        // Signify that the voice is active
      }                                     // ~~~~~~~~~ ResetState ~~~~~~~~~ //

      bool Process(
        AudioBuffer<double>& io,              // Output audio buffer
        size_t nf) noexcept                   // Number of frames to process
      {                                       // ~~~~~~~~ Process ~~~~~~~~~~~~ //
         auto* outL=io.Channel(0);            // Pointer to the left channel
         auto* outR=io.Channel(1);            // Pointer to the right channel
         for (size_t i=0;i<nf;++i)            // For each frame...
         {
            string->Propagate(1);             // Propagate the string for one sample
            double raw=static_cast<double>(string->Output()*gain); // Get the output of the string
            double s=tone.ProcessSample(raw); // Process the output through the tone shaping filter
            outL[i]+=s;                       // Add the processed sample to the left channel
            outR[i]+=s;                       // Add the processed sample to the right channel
            ++age;                            // Increment the age of the voice
         }                                    // Done writing and processing samples.
         active=age<static_cast<size_t>(fs*5.0);// 5 seconds hard kill.
         return active;                       // Return true if the voice is still active
      }                                       // ~~~~~~~~~ Process ~~~~~~~~~~~~ //

      void PitchBend(float semitones,size_t /*nframes*/) noexcept
      {                                       // ~~~~~~~~~ PitchBend ~~~~~~~~~~ //
        string->SetPitchBendCents(semitones*100); // Set the pitch bend in cents
      }                                       // ~~~~~~~~~ PitchBend ~~~~~~~~~~ //
      // Setters and getters
      inline void SetGain(double g) noexcept
      {
        gain=g;                              // Set the gain factor
      }
      inline double GetGain(void) const noexcept
      {
        return gain;                         // Get the gain factor
      }
      inline void SetLoopCutoff(double fc) noexcept
      {
        if (fc>=0.5*fs) return;
        loopfc=fc;                           // Set the loop loss filter cutoff frequency
        string->Prepare(fs,f0,loopfc,loopq,stage,D,mut,muf,order,pos); // Prepare the string element
      }
      inline double GetLoopCutoff(void) const noexcept
      {
          return loopfc;                       // Get the loop loss filter cutoff frequency
      }
      inline void SetLoopQ(double q) noexcept
      {
        if (q<=0.0) return;
        loopq=q;                            // Set the loop loss filter Q factor
        string->Prepare(fs,f0,loopfc,loopq,stage,D,mut,muf,order,pos); // Prepare the string element
      }
      inline double GetLoopQ(void) const noexcept
      {
        return loopq;                       // Get the loop loss filter Q factor
      }
      inline void SetToneCutoff(double fc) noexcept
      {
        if (fc>=0.5*fs) return;
        tonefc=fc;                          // Set the tone filter cutoff frequency
        Assemble();                         // Re-assemble the tone shaping filter
      }
      inline double GetToneCutoff(void) const noexcept
      {
        return tonefc;                      // Get the tone filter cutoff frequency
      }
      inline void SetToneQ(double q) noexcept
      {
        toneq=q;                            // Set the tone filter Q factor
        Assemble();                         // Re-assemble the tone shaping filter
      }
      inline double GetToneQ(void) const noexcept
      {
        return toneq;                       // Get the tone filter Q factor
      }
      inline void SetToneType(ToneType t) noexcept
      {
        tt=t;                               // Set the tone type
        Assemble();                         // Re-assemble the tone shaping filter
      }
      inline ToneType GetToneType(void) const noexcept
      {
        return tt;                          // Get the tone type
      }
      inline void SetIntegerDelay(double d) noexcept
      {
        D=d;                                 // Set the integer delay in samples
        string->SetIntegerDelay(D);         // Set the integer delay in the string element
      }
      inline double GetIntegerDelay(void) const noexcept
      {
        return D;                           // Get the integer delay in samples
      }
      inline void SetFractionalDelay(double m) noexcept
      {
        mut=m;                              // Set the Thiran fractional delay in samples
        string->SetFractionalDelay(mut);    // Set the Thiran fractional delay in the string element
      }
      inline double GetFractionalDelay(void) const noexcept
      {
        return mut;                         // Get the Thiran fractional delay in samples
      }
      inline void SetMuFarrow(double m1) noexcept
      {
        muf=m1;                             // Set the Farrow fractional delay in samples
        string->SetMuFarrow(muf);           // Set the Farrow fractional delay in the string element
      }
      inline void SetPositionEQ(double p) noexcept
      {
        pos=p;
        string->SetPositionEQ(pos);
      }
      inline double GetPositionEQ(void) const noexcept
      {
        return pos;                         // Get the position in the string
      }
      inline void SetScale(double s) noexcept
      {
        scale=s;                            // Set the scale factor for the pluck
      }
      inline void SetNoisePower(double p) noexcept
      {
        string->SetNoisePower(p); // Set the noise power for the string element
      }
      inline double GetNoisePower(void) const noexcept
      {
        return string->GetNoisePower(); // Get the noise power for the string element
      }
      inline void SetDispersionCoefficient(double a) noexcept
      {
        alpha=a;                            // Set the dispersion coefficient
        string->SetDispersionCoefficient(alpha); // Set the dispersion coefficient in the string element
      }
      inline double GetDispersionCoefficient(void) const noexcept
      {
        return alpha;                       // Get the dispersion coefficient
      }
      inline void SetDispersionStages(size_t st) noexcept
      {
        stage=st;                            // Set the number of dispersion stages
        string->SetDispersionStages(stage);  // Set the number of dispersion stages in the string element
      }
      inline size_t GetDispersionStages(void) const noexcept
      {
        return stage;                       // Get the number of dispersion stages
      }
      inline void SetBridgeFilter(double fc, double q) noexcept 
      {
        if (fc>=0.5*fs||q<=0.0) return;
       string->SetBridgeFilter(fc, q); // Set the bridge filter cutoff frequency and Q factor
      }
      inline void SetNutFilter(double fc, double q) noexcept 
      {
       if (fc>=0.5*fs||q<=0.0) return;
       string->SetNutFilter(fc, q); // Set the nut filter cutoff frequency and Q factor
      }
      inline void SetOrder(size_t o) noexcept
      {
        if (o<1) return;
        order=o;                            // Set the order of the Farrow deinterpolator
        string->SetOrder(order);            // Set the order in the string element
      }
      inline size_t GetOrder(void) const noexcept
      {
        return order;                       // Get the order of the Farrow deinterpolator
      }
      inline size_t AgeSamples(void) const noexcept
      {
        return age;                         // Get the age of the voice in samples
      }
    private:
      // ----------------------------------- //
      // Assemble the tone shaping filter
      // ----------------------------------- //
      inline void Assemble(void) noexcept
      {                                      // ~~~~~~~~~ Assemble ~~~~~~~~~ //
        switch(tt)                           // Dispatch filter according to desired tone type.
        {
            case ToneType::Lowpass: tone=ff.ButterworthLP(fs,tonefc,toneq); break; // Lowpass filter
            case ToneType::Highpass: tone=ff.ButterworthHP(fs,tonefc,toneq); break; // Highpass filter
            case ToneType::Bandpass: tone=ff.BandPass(fs,tonefc,toneq); break; // Bandpass filter
            case ToneType::Shelflow: tone=ff.LowShelf(fs,tonefc,-3.f); break; // Low shelf filter
            case ToneType::Shelfhigh: tone=ff.HighShelf(fs,tonefc,-3.f); break; // High shelf filter
            default: tone=ff.ButterworthLP(fs,tonefc,toneq); break; // Default to lowpass filter            
         }                                   // Done dispatching tone shaping filter
      }                                      // ~~~~~~~~~ Assemble ~~~~~~~~~ //
    private:
      // Members
      int cmidi{0};                          // Current MIDI note id.      
      double fs{48000.0};                    // Sample rate
      double f0{440.0};                      // Fundamental frequency  
      double loopfc{1000.0};                 // Loop loss filter cutoff frequency
      double loopq{0.7071};                  // Loop loss filter Q factor
      double tonefc{22000.0};                // Tone filter cutoff frequency
      double toneq{0.7071};                  // Tone filter Q factor
    
      double gain{1.0};                      // Gain factor
      double vel{0.0};                       // Velocity of the pluck
      double scale{0.4f};                    // Scale factor for the pluck
      double amp{0.0};                       // Amplitude of the pluck

      double D{0.0};                         // Integer delay in samples
      double mut{0.0};                       // Thiran fractional delay in samples
      double muf{0.0};                       // Farrow fractional delay in samples
      double alpha{0.18};                     // Dispersion coefficient
      size_t age{0};                         // Age of the voice
      size_t bs{256};                        // Block size of the chunk to process
      size_t ch{2};                          // Number of channels
      size_t samples{bs*ch};                 // Number of samples to process
      size_t stage{2};                       // Number of dispersion stages
      size_t order{K};                       // Order of the of the filters
      bool active{false};                     // True if the voice is active
      double pos{0.15};                      // Position in the string
      // Our objects
      std::unique_ptr<String> string;         // Pointer to the string element
      FilterFactory<double> ff;               // Filter factory for easy filter creation
      BiQuad<double> tone;                    // Tone shaper
      ToneType tt{ToneType::Lowpass};         // Tone type
  };
}