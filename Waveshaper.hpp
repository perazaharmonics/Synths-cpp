/*
* * 
* * Filename: Waveshaper.hpp
* *
* * Description:
* *   A Wave Shaper is a signal processing element that applies a non-linear transformation to an input signal.
* *   It can be used to create various waveforms, distortions, and effects.
* *   The WaveShaper class provides a static method to shape a sine signal based on the specified waveform ID.
* *   
* *
* *
* * Author:
* * JEP J. Enrique Peraza
*/
#pragma once
#include <vector>
#include <tuple>
#include <cmath>
#include <iostream>

namespace sig::wg
{
namespace detail{template<typename T>constexpr T TWOPI(void){return static_cast<T>(6.2831853071795864769);}}
/* ------------------------------------------------------------ */
/*  Shared helper: wave shaper &&  vibrato LFO                  */
/* ------------------------------------------------------------ */
template<typename T>
struct WaveShaper
{
  static inline int Hash(int x) noexcept
  {
    x^=x>>13;                         // XOR with right shift
    x*=0x85ebca6b;                    // Multiply by a large prime
    x^=x>>16;                         // XOR with another right shift
    return x;                         // Return the hashed value
  }                                   // ~~~~~~~~~~ Hash ~~~~~~~~~~ //
  static inline T Grad(
    int h,                             // Integer part
    T x) noexcept                       // fractional part
  {                                     // ~~~~~~~~~~~ Grad ~~~~~~~~~~~~ //
    return ((h&1)?x:-x);                // 1-D gradient function
  }                                     // ~~~~~~~~~~~ Grad ~~~~~~~~~~~~ // 
  static inline T Perlin(T x)           // Generate Perlin noise
  {                                     // using arbitrary point x in space.
    int xi=static_cast<int>(std::floor(x))&255; // Get the integer part of x
    T xf=x-std::floor(x);               // Get the fractional part of x
    T u=xf*xf*xf*(xf*(xf*6-15)+10);     // Smoothstep interpolation
    int h0=Hash(xi);                    // Hash the integer part
    int h1=Hash(xi+1);                  // Hash the next integer part
    T g0=Grad(h0,xf);                   // Get the gradient for the first
    T g1=Grad(h1,xf-1);                 // Get the gradient for
    return ((1-u)*g0+u*g1);             // Interpolate between the two gradients
  }                                     // ~~~~~~~~~~ Perlin ~~~~~~~~~~ //
  static T Shape(
    T s,                               // A sine signal 
    T ph,                              // Phase value
    int id,                            // Waveform ID
    std::mt19937& rng) noexcept        // Random number generator        
  {
    using std::abs; using std::asin; using std::exp; using std::log;
    switch (id&0xF)                     // Switch according to the wave ID.
    {
      case 1:return (2.0/M_PI)*asin(s); // Triangle wave.
      case 2:return (s>=0?)1:-1;        // Square wave.
      case 3:return ph/M_PI-1;          // Sawtooth wave.
      case 4: return 1-ph/M_PI;         // Inverted sawtooth wave.
      case 5: return exp(s)-1.718281828;// Exponential wave.
      case 6: return log(abs(s)+1e-3);  // Logarithmic wave.
      case 7: return s*s*(s>0?1:-1);    // Signed sine squared wave.
      case 8:return abs(s);             // Absolute value wave.
      case 9: return ph<M_PI?1:-0.5;    // Half sine wave.
      case 10: return ph<M_PI?ph/M_PI:-1;// Sharkfin wave.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Sample & Hold waveform
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      case 11:
      {
          static thread_local T h=T(0); // Thread-local variable for noise
          static thread_local int c=0;  // Thread-local counter
          if (++c>128)                  // Counter exceeds half the block size?
          {                             // Yes, generate new noise
            std::uniform_real_distribution<T> dist(-1.0, 1.0);
            h=dist(rng);                // Generate a new noise distribution
            c=0;                        // Reset the counter
          }                             // Done generating new noise
          return h;                     // Return the current noise value
      }                                 // End of S&H waveform
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Smooth noise.                  //
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      case 12:                          // Smooth noise
      {                                 // 
        static thread_local T y=0;      // Thread-local variable for noise
        std::uniform_real_distribution<T> dist(-1.0, 1.0);
        y+=(dist(rng)-y)*0.5;           // Smooth the noise value
        return y;                       // Return the current noise value
      }                                 // End of smooth noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Pink noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      case 13:
      {
        static thread_local T b0{},b1{},b2{};// Thread-local variables for pink noise
        std::uniform_real_distribution<T> dist(-1.0, 1.0);
        T w=dist(rng);                  // Generate a new random value
        b0=0.997*b0+0.029591*w;         // Update first filter coeff
        b1=0.985*b1+0.032534*w;         // Update second filter coeff
        b2=0.95*b2+0.048056*w;          // Update third filter coeff
        return (b0+b1+b2+w*0.05)*0.25;  // Return the pink noise.
      }                                 // Done with pink noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Worley noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      case 14:  
      {       
        static thread_local T v=0;      // Thread-local variable for Worley noise
        static thread_local int c=0;    // Thread-local counter
        if (++c>256)                    // Counter exceeds block size?
        {                               // Yes, generate new noise
          std::uniform_real_distribution<T> dist(0.0, 1.0);
          v=dist(rng);                  // Generate a new noise distribution
          c=0;                          // Reset the counter
        }                               // Done generating new noise
        return v;                       // Return the current noise value
      }                                 // Done with Worley noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Perlin noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      case 15:
      {                                 // Perlin noise
        static thread_local T t=0;      // Thread-local variable for Perlin noise
        t+=.25;                         // Increment the Perlin noise time.
        return Perlin(t)*0.8;           // Scale the Perlin noise and return it.   
      }                                 // Done with Perlin noise.
      default: return s;                // Default case,return sine wave.
    }                                   // Done with dispatcher
  }                                     // ~~~~~~~~~~ Shape ~~~~~~~~~~ //
};
}
