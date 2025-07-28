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
// ------------------------------------------------------------ //
// Wave Shaping dispatcher.                 
// ------------------------------------------------------------ //
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
  static T Simplex1D(T x) noexcept 
  {                                     // ~~~~~~~~~ Simplex1D ~~~~~~~~~ //
    int i=static_cast<int>(std::floor(x));// Get the integer part of x
    T f=x-i;                            // Get the fractional part of x
    T g0=Grad(Hash(i+1),f-1);           // Get the gradient for the first
    T t0=0.5-f*f;                       // Compute the first interpolation factor
    t0=t0<0?0:t0*t0*t0*t0;              // Smoothstep interpolation
    t1=0.5-(f-1)*(f-1);                 // Compute the second interpolation factor
    t1=t1<0?0:t1*t1*t1*t1;              // Smoothstep interpolation
    return T(8)*(t0*g0+t1*g1);          // Return the final value
  }                                     // ~~~~~~~~~ Simplex1D ~~~~~~~~~ //
  static T FBM(T x) noexcept            // Fractal Brownian Motion Noise
  {                                     // ~~~~~~~~~ FBM ~~~~~~~~~~~~~~~ //
    T a=1;                              // Initial amplitude
    T s=0;                              // Initial sum
    for (int o=0;o<5;++o)               // For the number of octaves
    {                                   // Generate the fractal brownian noise
      s+=a*Perlin(x);                   // Add the Perlin noise value
      x*=2;                             // Double the frequency
      a*=0.5;                           // Halve the amplitude
    }                                   // Done generating FBM
    return 0.7*s;                       // Scale the final value
  }                                     // ~~~~~~~~~ FBM ~~~~~~~~~~~~~~~ //
  static T WaveletNoise(T x) noexcept
  {                                     // ~~~~~~~~~ WaveletNoise ~~~~~~ //
    return (Perlin(x)-0.5*Perlin(2*x))*0.67;// Generate wavelet noise
  }                                     // ~~~~~~~~~ WaveletNoise ~~~~~~ //
  static T Shape(
    T s,                               // A sine signal 
    T ph,                              // Phase value
    int id,                            // Waveform ID
    std::mt19937& rng) noexcept        // Random number generator        
  {
    using std::abs; using std::asin; using std::exp; using std::log;
    switch (id%20)                      // Switch according to the wave ID.
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
          if (++c>128)                  // Counter exceeds the grid size?
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
        if (++c>256)                    // Counter exceeds grid size?
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
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Fractal Brownian Motion Noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      case 16:
      {
        static thread_local T t=0;      // Thread-local variable for FBM noise
        t+=.25;                         // Increment the FBM noise time.
        return FBM(t);                  // Return the scaled FBM noise.
      }                                 // Done with FBM noise.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Simplex1D noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      case 17:
      {
        static thread_local T t=0;      // Thread-local variable for Simplex1D noise
        t+=.25;                         // Increment the Simplex1D noise time.
        return Simplex1D(t);            // Return the scaled Simplex1D noise.
      }                                 // Done with Simplex1D noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Wavelet noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      case 18:
      {
        static thread_local T t=0;      // Thread-local variable for Wavelet noise
        t+=.25;                         // Increment the Wavelet noise time.
        return Wavelet(t);              // Return the scaled Wavelet noise.
      }                                 // Done with Wavelet noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // OpenSimplex1D Proxy noise
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      case 19:
      {
        static thread_local T t=0;      // Thread-local variable for OpenSimplex1D noise
        t+=.25;                         // Increment the OpenSimplex1D noise time.
        return OpenSimplex1D(t*1.5);    // Return the scaled OpenSimplex1D noise.
      }
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Default case
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      default: return s;                // Default case,return sine wave.
    }                                   // Done with dispatcher
  }                                     // ~~~~~~~~~~ Shape ~~~~~~~~~~ //
};
} // namespace sig::wg
