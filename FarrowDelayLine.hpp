 /*
 * *
 * * Filename: FarrowDelayLine.hpp
 * *
 * * Description:
 * * A time-varying fractional delay using Farrow's method for Lagrange interpolation.
 * *  instead of using the Direct Form FIR filter, we use a Farrow structure
 * *  to compute the fractional delay. A Farrow filter rewrites any polynomial
 * *  based fractional delay FIR as a sum of fixed sub-filter whore outputs are
 * *  combined by powers of mu. It's nice and very stable for audio-rate modulation
 * *  of delay lines, like pitch bends. It's less compuationally consumptive too.
 * *     y[n, mu] = sum_{m=0}^M (mu^m * v_m[n]), with v_m[n] = sum_{k=0}^K (c_{m,k} * x[n - D - k])
 * * 
 * * Where:
 * *   M = polynoimial order (for an N-th order Lagrance Fractional Delay M=N).
 * *   c_{m,k} are the fixed FIR taps (constant at run-time).
 * *   All time-varying cost is in the Horner-style evaluation of sum_{m=0}^M (mu^m * v_m[n]).
 * * Author:
 * *  JEP J.Enrique Peraza
 * *
  */
#pragma once
#include <array>
#include <cstddef>
#include <algorithm>
#include <cmath>
#include "DelayLine.hpp"

namespace sig::wg
{
  // ---------------------------------- //
  // Helper: Compile-time 2-D array for coefficients matrix (max order +1)^2
  // ---------------------------------- //
  template<typename T, size_t MaxOrder>
  using CoeffMatrix=std::array<std::array<T,MaxOrder+1>,MaxOrder+1>;

  // ---------------------------------- //
  // Polynomial helpers - naive because orders are tine (<=5)
  // ---------------------------------- //
  namespace detail
  {
    // -------------------------------- //
    // Expand product (x-m) for all m in the provided list, normalize by
    // the denominator.
    // -------------------------------- //
    template <typename T, size_t Max>
    void MultiplyMonomial(              // Multiply a monomial (x-m) for all m in the provided list
      std::array<T,Max>& p,               // Polynomial coefficients array
      size_t order,                // Order of the polynomial
      T root,                           // The root to multiply by
      T denom) noexcept                 // The denominator to normalize by
    {                                   // ----------- MultiplyMonomial ----------------- //
        // ---------------------------- //
        // Poly holds coefficients up to current degree (inclusive)
        // newpoly=poly*(x-root)/denominator
        // ---------------------------- //
        for (size_t i=std::min<std::size_t>(order+1,Max-1);i-->0;)  // For each coefficient in the polynomial.
        {                               // Compute new polynomial coefficients.
            p[i+1]+=p[i]/denom;         // Compute the new coefficient.
            p[i]*=(-root)/denom;        // Update the next coefficient.
        }                               // Done with the new polynomial coefficients.
    }                                   // ----------- MultiplyMonomial ----------------- //
    // -------------------------------- //
    // Build Lagrance basis polynomial L-k(x) of degree N (0<=k<=N)
    // -------------------------------- //
    template <typename T,size_t Max>
    std::array<T,Max> BuildLagrangeCoeffs(
      size_t k,                         // The index of the Lagrange polynomial to build
      size_t N) noexcept                // The order of the Lagrange polynomial
    {                                   // ----------- BuildLagrangeCoeffs ------------- //
      std::array<T,Max> coefs{}; // Coefficients of the Lagrange polynomial
      coefs[0]=T(1);                    // Initialize the first coefficient to 1
      size_t deg=0;                     // Degree of the polynomial 
      for (size_t m=0;m<=N;++m)         // For the order of the polynomial
      {                                 // Compute the coefficients of the Lagrange polynomial
        if (m==k)                       // Is the current index equal to K?
          continue;                     // Yes not this one.
        MultiplyMonomial(coefs,deg,static_cast<T>(m),static_cast<T>(k)-static_cast<T>(m));// Multiply the monomial (x-m) for all m in the provided list
        ++deg;                          // Increment the degree of the polynomial
      }                                 // Done with the coefficients of the Lagrange polynomial
      return coefs;                     // Return the coefficients of the Lagrange polynomial   
    }                                   // ----------- BuildLagrangeCoeffs ------------- //
  }                                     // namespace detail
  // ---------------------------------- //
  // Farrow-structure Lagrange Fractional-Delay Filter (order <= MaxOrder)
  // ---------------------------------- //
  template<typename T=float,size_t MaxOrder=5>
  class FarrowFDFilter
  {
    public:
      constexpr static size_t MAXORDER=MaxOrder;
      ~FarrowFDFilter(void) noexcept = default; // Default destructor
      void SetOrder(
        size_t N) noexcept              // Set the order of the Lagraange Interpolator
      {                                 // ----------- SetOrder ----------------- //                 
        order=std::min<size_t>(N,MaxOrder); 
        BuildCoefficients();           // Build the coefficients for the Lagrange filter
      }                                // ----------- SetOrder ----------------- //
    void SetMu(                      // Set the fractional part of the delay line
       T mu ) noexcept                // The fractional part of the delay 
    {                                  // ----------- SetMu ----------------- //
      this->mu=mu-std::floor(mu);      // Set the fractional part of the delay
    }
    template<size_t Len>             // Length of the delay line.
    bool Process(
      const DelayLine<T,Len>& dl,   // The delay line to process
      size_t D,                     // The fractional delay in samples
      T* const y) noexcept      // Output buffer.
    {                               // ---------- Process -------------- //
      size_t maxD=Len-1-order;   // Maximum delay length.
      if (D>maxD) D=maxD;              // Clamp the delay to the
      // -------------------------- //
      // v_m=S_k (C[m][k]*x[n-D-k]))
      // -------------------------- //
      std::array<T,MaxOrder+1> v{}; // Zero the output buffer
      for (size_t k=0;k<=order;++k) // For each coefficient in the Lagrange polynomial
      {
        const T xk=dl.Peek(D+k);    // Get the sample at index D+k from the delay line
        for (size_t m=0;m<=order;++m)
          v[m]+=coeff[m][k]*xk;    // Compute the output sample by summing the products of the coefficients and the corresponding samples from the delay line
      }                            // Done with the output samples.
      // ------------------------- //
      // Perform Horner's evalueation of coeffiecient expansion
      // y=(((v_N)*mu+n_{N-1})*mu+ ... +v_0)*mu+v_0
      // ------------------------- //
      *y=v[order];                 // Initialize the output sample with the last coefficient
      for (size_t m=order;m-->0;)  // For each coefficient in the Lagrange polynomial
        *y=(*y)*mu+v[m];           // Compute the output sample by summing the products of the coefficients and the corresponding samples from the delay line
      return true;                 // Return true to indicate success
    }                              // ---------- Process -------------- //
    private:
    size_t order{3};               // Order of the Lagrange filter, default is 3.
    T mu=T(0);                     // Fractional part of the delay, default is 0.
    CoeffMatrix<T,MaxOrder> coeff{}; // Coefficients for the Lagrange filter, size is order+1.
    void BuildCoefficients(void) noexcept // Build the coefficients for the Lagrange filter
    {                                   // ----------- BuildCoefficients ----------------- //
      // ------------------------------ //
      //  coeff[m][k]=polynomial coefficients of mu^m for tap k
      // ------------------------------ //
      for (size_t k=0;k<=order;++k)     // For each coefficient in the interpolator
      {
        auto c=detail::BuildLagrangeCoeffs<T,MaxOrder+1>(k,order);
        for (size_t m=0;m<=order;++m)   // For the order of the filter.
          coeff[m][k]=c[m];             // Store the coefficients in the matrix
      }                                 // Done with the coefficients.                   
    }                                   // ----------- BuildCoefficients ----------------- //
  };                                    // class FarrowFDFilter
  // ---------------------------------- //
  // Variable fractional delay line using Farrow's approach to Lagrange interpolation
  // ---------------------------------- //
  template<typename T=float,
    size_t MaxLen=1024,
    size_t MaxOrder=5>
  
  class FarrowDelayLine
  {
    public:
      FarrowDelayLine(void) noexcept
      {
        dl=new sig::DelayLine<T,MaxLen>{}; // Create a new delay line with the specified maximum length.
        fd=new FarrowFDFilter<T,MaxOrder>{}; // Create a new Farrow filter with the specified maximum order.
        fd->SetOrder(order);               // Set the order of the Lagrange interpolator.
        delay=T(1.3);                        // Initialize the delay to 0.
        const float m=delay-std::floor(delay);
        this->mu=m;
        fd->SetMu(m);                   // Set the fractional part of the delay line to 0.
        
        incr=T(0.01);                         // Initialize the increment to 0.
        targ=std::floor(delay+m);                         // Initialize the target delay to 0.
      }
      ~FarrowDelayLine(void)
      {
        delete fd;
        fd=nullptr;
        delete dl;
        dl=nullptr;
      }
      T ReadFrac(
        T d) const noexcept                   // The fractional delay
      {
        return dl->ReadFrac(d);         // Read the sample from the delay line using linear interpolation.
      }
      inline void Write(const T& sample) noexcept
      {
        dl->Write(sample);              // Write the sample to the delay line.
        haswritten=true; // Mark that we have written a sample.
      }
      void Prepare(                     // Prepare the delay line for processing 
        T d)                            // The integer delay
      noexcept
      {                                 // ----------- Prepare ----------------- //
        if (d<=T(0.0)) d=T(0.0); // Ensure the delay is at least 1 to avoid underflow.
        else if (d>=T(MaxLen-1)) d=T(MaxLen-1); // Ensure the delay is at most 1 to avoid overflow.
        assert(d>=T(0.0)&&d<=T(MaxLen-1));   // Ensure the delay is within bounds.
        order=std::min<size_t>(order,MaxOrder); // Set the order of the interpolator.
        fd->SetOrder(order);             // Set the order
        dl->Clear();
        fd->SetOrder(order);              // Set the order of the Lagrange interpolator.
        this->delay=d;                // Set the desired delay.
        const float m=delay-std::floor(delay);        
        fd->SetMu(m);                  // Set the fractional part of the delay line to 0.
        targ=std::floor(delay+mu);          // Set the target delay to the integer part of the delay plus the fractional part.
      }                                 // ----------- Prepare ----------------- //
      void Propagate(size_t n) noexcept
      {                                 //% Circulate n sampples through delay line.
        for (size_t i=0;i<n;++i)        // For each sample in the block
          Tick();                       // Call Tick to update the delay line and interpolator state.
      }
      void Clear(void) noexcept
      {                                 // ----------- Clear ----------------- //
        dl->Clear();                     // Clear the delay line buffer.
        fd->SetOrder(order);             // Set the order of the Lagrange interpolator.
        fd->SetMu(T(0));                 // Set the fractional part of the delay line to 0.
        haswritten=false; // Mark that we have not written any samples.
      }                                 // ----------- Clear ----------------- //

    // Read interpolated output
    T Read(void) noexcept
    {                                   // ----------- Read ------------------ //
      
      size_t D=static_cast<size_t>(std::floor(delay)); // Get the integer part of the delay.
      T mu=delay-D; // Get the fractional part of the delay.
      if (mu==T(0))                    // No fractional delay?
        return dl->Peek(D+1);
      if (delay<T(1))                  // Fractional delay?
        return dl->ReadFrac(delay); // Read the sample from the delay line using linear interpolation.
      // Otherwise, run the Farrow filter
      T y{};                            // Output sample
      if (!fd->Process(*dl,D,&y))       // Process the delay line with the Farrow filter
        y=T(0);                         // Return 0 if processing fails.
      return y;                         // Return the output sample.
    }                                   // ----------- Read ------------------ //
    // Smooth glide: move to a new delay in k samples.
    void RampTo(
      T targ,                           // Target delay to ramp to
      size_t k) noexcept                // Number of samples to ramp to the target
    {                                   // ---------- RampTo ----------------- //
      if (targ<T(0.0)) targ=T(0.0); // Ensure the target delay is at least 1 to avoid underflow.
      else if (targ>=T(MaxLen-1)) targ=T(MaxLen-1); // Ensure the target delay is at most 1 to avoid overflow.
      assert(targ>=T(0.5)&&targ<=T(MaxLen-1));
      if (delay<T(0.0)) delay=T(0.0); // Ensure the current delay is at least 1 to avoid underflow.
      else if (delay>=T(MaxLen-1)) delay=T(MaxLen-1); // Ensure the current delay is at most 1 to avoid overflow.
      assert(delay>=T(0.5)&&delay<=T(MaxLen-1)); //
      incr=(targ-delay)/static_cast<T>(k); // Compute the increment for the ramp.
      this->targ=targ;                  // Set the target delay.
    }                                   // ---------- RampTo ----------------- //
    // Call once per sample to progress an active glide.
    void Tick(void) noexcept
    {                                   // ---------- Tick ------------------- //
      // Something happened... update myself.
      if((incr>0&&delay<targ)||(incr<0&&delay>targ))
      {                                 // We updated the incr or delay greater than target?
        delay+=incr;                    // We increased delay by this much.
        // If we are ramping push delay below 1, clamp targ as well
      if (delay<T(0.0)) delay=T(0.0); // Ensure the delay is at least 1 to avoid underflow.
      else if (delay>=T(MaxLen-1)) delay=T(MaxLen-1);
      assert(delay>=T(0.0)&&delay<=T(MaxLen-1)); // Ensure the delay is within bounds.
      if (targ<T(0.0)) targ=T(0.0); // Ensure the target delay is at least 1 to avoid underflow.
      else if (targ>=T(MaxLen-1)) targ=T(MaxLen-1); // Ensure the target delay is at most 1 to avoid overflow.
      assert(targ>=T(0.0)&&targ<=T(MaxLen-1)); // Ensure the target delay
      const float m=delay-std::floor(delay); // Get the fractional part of the delay.
      fd->SetMu(m); // Update the interpolator with the delay state.
      }                                 // Done updating out state.
    }                                   // ---------- Tick ------------------- //
    // -------------------------------- //
    // Set exact (possibly non-integer) delay 'd' with Lagrange Order 'N'.
    // Exists in the range of 0<= d < MaxLen-1.
    // -------------------------------- //
    void SetDelay(
      T d,                              // Fractioal delay to ramp to.
      size_t N=3) noexcept              // Default filter order
    {                                   // ----------- SetDelay ---------------- //
      if (d<=T(0.0)) d=T(0.0); // Ensure the delay is at least 1 to avoid underflow.
      else if (d>=T(MaxLen-1)) d=T(MaxLen-1); // Ensure the delay is at most 1 to avoid overflow.
      order=std::min<std::size_t>(N,MaxOrder);// Set the order of the interpolator.
      fd->SetOrder(order);               // Set the order of the Lagrange interpolator.
      delay=d;                          // Set the desired delay.
      const float m=delay-std::floor(delay); // Get the fractional part of the delay.
      this->mu=m;                      // Set the fractional part of the delay line.
      fd->SetMu(m); // Set the fractional part of the delay.
    }                                   // ----------- SetDelay ---------------- //
    void SetMu(
      T mu) noexcept                    // Set the fractional part of the delay line
    {                                   // ---------- SetMu ----------------- //
      if (mu<=-T(MaxLen-1)) // Ensure the mu is at least 0.
        mu=-(T(MaxLen-1)); // If not, set it to the fractional part of the delay.
      else if (mu>=T(MaxLen-1))             // Ensure the mu is at most 1.
        mu=T(MaxLen-1);                     // If not, set it to 1.
      this->targ+=std::floor(mu+delay);
      this->mu=mu;
      fd->SetMu(mu);                    // Set the fractional part of the delay line.
    }          
    void SetIncrement(
      T incr) noexcept                  // Set the increment for the ramp
    {                                   // ---------- SetIncrement ----------------- //
      this->incr=incr;                  // Set the increment for the ramp.
    }                                   // ---------- SetIncrement ----------------- //

    private:
      sig::DelayLine<T,MaxLen>* dl{nullptr}; // Delay line buffer
      FarrowFDFilter<T,MaxOrder>* fd{nullptr}; // Farrow filter for
      size_t order{3};                // Order of the Lagrange filter, default is 3.
      T mu{0.0f};
      T delay{1.0f};                     // Current delay in samples.
      T incr{0.01};                      // Increment for the ramp.
      T targ{0.0f};  // The integer delay                      // Target delay to ramp to.
      bool haswritten{false};
  };                                  // class FarrowDelay
}                                     // namespace sig::wg
