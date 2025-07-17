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
  template<typename T=float,
  size_t MaxLen=1024,
  size_t MaxOrder=5>
  class FarrowInterpolator
  {
    public:
      constexpr static size_t MAXORDER=MaxOrder;
      ~FarrowInterpolator(void) noexcept = default; // Default destructor
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
    
    bool Process(
      const DelayLine<T,MaxLen>& dl,   // The delay line to process
      size_t D,                     // The fractional delay in samples
      T* const y) noexcept      // Output buffer.
    {                               // ---------- Process -------------- //
      if (D<order || D > MaxLen - 1) return false;
      size_t maxD=MaxLen-1-order;   // Maximum delay length.
      if (D>maxD) D=maxD;              // Clamp the delay to the
      // -------------------------- //
      // v_m=S_k (C[m][k]*x[n-D-k]))
      // -------------------------- //
      std::array<T,MaxOrder+1> v{}; // Zero the output buffer
      for (size_t k=0;k<=order;++k) // For each coefficient in the Lagrange polynomial
      {
        const T xk=dl.Peek(D-k);    // Get the sample at index D+k from the delay line
        for (size_t m=0;m<=order;++m)
          v[m]+=coeff[m][k]*xk;    // Compute the output sample by summing the products of the coefficients and the corresponding samples from the delay line
      }                            // Done with the output samples.
      // ------------------------- //
      // Perform Horner's evalueation of coeffiecient expansion
      // y=(((v_N)*mu+n_{N-1})*mu+ ... +v_0)*mu+v_0
      // ------------------------- //
      *y=0;                 // Initialize the output sample with the last coefficient
      for (int i=static_cast<int>(order);i>=0;--i)  // For each coefficient in the Lagrange polynomial
        *y=(*y)*mu+v[i];           // Compute the output sample by summing the products of the coefficients and the corresponding samples from the delay line
      return true;                 // Return true to indicate success
    }                              // ---------- Process -------------- //
    inline const std::array<T, MaxOrder+1>& operator[](size_t i) const noexcept // ----------- GetCoefficients ----------------- //
    {
      assert(i < MaxOrder+1);
      return coeff[i];
    }
    inline std::array<T, MaxOrder+1>& operator[](size_t i) noexcept // ----------- GetCoefficients ----------------- //
    {
      assert(i < MaxOrder+1);
      return coeff[i];
    }
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
  };                                    // class FarrowInterpolator
  // ---------------------------------- //
  // Farrow Deinterpolator to estimate the fractional delay of where to
  // write the sample in the delay line. As opposed to the interpolator, that determins
  // what the output sample iss, the deinterpolator determines where the input sample
  // should be written in the delay line.
  // ---------------------------------- //
  template<typename T=float,
    size_t MaxLen=1024,
    size_t MaxOrder=5>
  class FarrowDeinterpolator
  {
    public:
      FarrowDeinterpolator(void) noexcept
      {
        BuildCoefficients();           // Build the coefficients for the Lagrange filter
      }
      ~FarrowDeinterpolator(void) noexcept = default; // Default destructor
      /// Note that the order MUST be the same as in the interpolator.
      void SetOrder(size_t N) noexcept // Set the deinterpolator order
      {
        order=std::min<size_t>(N,MaxOrder);// Clamp the order to the maximum order
        BuildCoefficients();           // Build the coefficients for the Lagrange filter
      }                                // ----------- SetOrder ----------------- //
      void SetMu(T m) noexcept        // Set the fract part of the dela
      {
        mu=std::clamp(m-std::floor(m),T(0),std::nextafter(T(1),T(0)));
      }
      T GetMu(void) const noexcept    // Get the fractional part of the delay
      {
        return mu;                     // Return the fractional part of the delay
      }
      // Apply the inverse Farrow filter one write (distribute x into fractional slots)
      bool Process(
        T x,                          // The sample to write to the delay line
        DelayLine<T,MaxLen>* dl,      // The delay line to write to (mutable)
        size_t D) noexcept            // The total delay in samples
        {
          if (D<order || D>MaxLen-1) return false;
          std::array<T,MaxOrder+1> v{};   // Zero the output buffer
          // ------------------------- //
          // Calculate the coefficients using Horner's method fo mu^m
          // ------------------------- //
          for (size_t m=0;m<=order;++m)// Rows
          {                            
            for (size_t k=0;k<=order;++k)// Cols
            {
              if (k==0)
                v[m]=coeff[m][k]*x; // Initialize the first coefficient
              else
                v[m]+=coeff[m][k]*std::pow(mu,k); // Compute the coefficient for the current order and tap
            }
          }                              // Done with the coefficients.
          // ------------------------- //
          // Perform additive write to the delay line
          // ------------------------- //
          for (size_t k=0;k<=order;++k) // For each coefficient in
            dl->Write(D-k,dl.Peek(D-k)+v[k]);// Additiive write to the delay line
          return true;                 // Return true to indicate success
        }                                // ---------- Process -------------- //
    private:
      size_t order{3};                // Order MUST be the same as in the interpolator.
      T mu{0};                        // Fractional part of the delay, default is 0.
      CoeffMatrix<T,MaxOrder> coeff{}; // Coefficients for the Lagrange filter, size is order+1.
      void BuildCoefficients(void) noexcept
      {
        for (size_t k=0;k<=order;++k)
        {                              // For each coefficient in the filter...
          auto c=detail::BuildLagrangeCoeffs<T,MaxOrder+1>(k,order);
          for (size_t ,=0;m<=order;++m)
            coeff[m][k]=c[m];         // Store the coefficients in the matrix
        }
      }
  };                                  // class FarrowDeinterpolator
  
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
        fd=new FarrowInterpolator<T,MaxLen,MaxOrder>{}; // Create a new Farrow filter with the specified maximum order.
        fdi=new FarrowDeinterpolator<T,MaxLen,MaxOrder>{}; // Create a new Farrow deinterpolator with the specified maximum order.
        fd->SetOrder(order);               // Set the order of the Lagrange interpolator.
        fdi->SetOrder(order);              // Set the order of the Lagrange deinterpolator.
        mu=std::clamp(delay - std::floor(delay), T(0.0), std::nextafter(T(1.0), T(0.0)));
        fd->SetMu(mu);                   // Set the fractional part of the delay line to 0.
        fdi->SetMu(mu);                  // Set the fractional part of the delay line to 0.
        incr=T(0.01);                         // Initialize the increment to 0.
        targ=std::floor(delay+mu);                         // Initialize the target delay to 0.
      }
      ~FarrowDelayLine(void)
      {
        delete fd;
        fd=nullptr;
        delete dl;
        dl=nullptr;
        delete fdi;
        fdi=nullptr;
      }
      T ReadFrac(
        T d) noexcept                   // The fractional delay
      {
        const T mu=d-std::floor(d); // Get the fractional part of the delay.
        const size_t D=static_cast<size_t>(std::floor(d));
        const size_t i=dl->GetHead(); // Write head index
        // ----------------------------  //
        // Access samples around index i: s[i], s[i=1], s[i+2], ...
        // ----------------------------  //
        T acc=0;
        constexpr size_t Order=MaxOrder; // Order of the Lagrange interpolator
        for (size_t k=0; k<=Order;++k)  // For coeffs of this order
        {                               // circuilate
          const T coeff=(*fd)[k][0]+
            (*fd)[k][1]*mu+      // Fourth order polynomial interpolation
            (*fd)[k][2]*mu*mu+    //
            (*fd)[k][3]*mu*mu*mu;  // ************************************* Circulate 
          acc+=coeff*dl->Peek(i-D-k); // Accumulate the result by multiplying the coefficient with the sample at the specified index.
        }                              // Done with the coefficients.
        return acc;                   // Return the accumulated result.
      }
      // Read reverse assumes we reverse relative to current delay line.
      T ReadReverse(T d) noexcept       // Delay with int and frac part
      {                                 //%ReadReverse:
        const size_t Dbase=static_cast<size_t>(std::floor(d)); // Get the integer part of the delay.
        const size_t D=(Dbase>static_cast<size_t>(std::floor(d)))?Dbase-static_cast<size_t>(std::floor(d)):0;
        const T revm=mu+(d-std::floor(d));// Reverse direction so add mu.
        const size_t revD=(revm>=T(1.0))?D+1:D; // Get the reverse delay.
        const T fmu=(revm>=T(1.0))?revm-T(1.0):revm;// Final fractional part of the delay.
        fd->SetMu(fmu); // Set the fractional part of the delay in the Farrow filter.
        T y{};                          // Output sample
        if (!fd->Process(*dl,revD,&y))  // Process the delay
          return T(0);                  // Return 0 if processing fails.
        return y;                       // Return the output sample.
      }                                 //%ReadReverse.
      
      inline void Write(const T& sample) noexcept
      {
        dl->Write(sample);              // Write the sample to the delay line.
        haswritten=true; // Mark that we have written a sample.
      }
      inline void WriteFrac(
        T x,                            // The sample to write to the delay line
        T d) noexcept                   // Delay with int and fractional part
      {
        d=std::clamp(d,T(0),T(MaxLen-1)); // Clamp the delay to the range [0, MaxLen-1].
        const size_t D=static_cast<size_t>(std::floor(d));// Get int part of the delay.
        const T m=d-T(D);               // Get the fractional part of the delay.
        fdi->SetMu(m);                  // Set the fractional part of the delay in the deinterpolator.
        if (!fdi->Process(x,dl,D))           // Process delay line with the deinterpolator.
          return;
        haswritten=true; // Mark that we have written a sample.
      }
      void Prepare(                     // Prepare the delay line for processing
        T d)                            // The delay (may be fractional, whole, both or int)
      noexcept
      {                                 // ----------- Prepare ----------------- //
        delay=std::clamp(d,T(0.0),T(MaxLen-1)); // Clamp the delay to the range [0, MaxLen-1].
        // Extract fractional part of the delay line
        mu=delay-std::floor(delay);
        // Rebuild farrow filter
        fd->SetOrder(order);           // Set the order of the Lagrange interpolator.
        fdi->SetOrder(order);          // Set the order of the Lagrange deinterpolator.
        fd->SetMu(mu);                 // Set the fractional part of the delay line.
        fdi->SetMu(mu);                // Set the fractional part of the delay line.
        dl->Clear();                 // Clear the delay line buffer.
        // set a matching targ so RampTo() won't fire until called.
        incr=T(0.0);                  // Reset the increment to 0.
        targ=delay;
        haswritten=false; // Reset the written flag.
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
        delay=T(0.0);                   // Reset the delay to 0.
        mu=T(0.0);                       // Reset the fractional part of the delay line to 0.
        incr=T(0.0);                     // Reset the increment to 0.
        targ=0;                          // Reset the target delay to 0.
        haswritten=false; // Mark that we have not written any samples.
      }                                 // ----------- Clear ----------------- //

    // Read interpolated output
    T Read(void) noexcept
    {                                   // ----------- Read ------------------ //
      
      size_t D=static_cast<size_t>(std::floor(delay)); // Get the integer part of the delay.
      T f=delay-D; // Get the fractional part of the delay.
      f=std::clamp(f,T(0.0),std::nextafter(T(1.0),T(0.0))); // [0,1) range
      if (D<order)
        return dl->Peek(0);   // Fallback most recent sample.
      if (std::abs(f)<=std::numeric_limits<T>::epsilon())                    // No fractional delay?
        return dl->Peek(D);
      // Otherwise, run the Farrow filter
      T y{};                            // Output sample
      if (!fd->Process(*dl,D,&y))       // Process the delay line with the Farrow filter
        y=T(0);                         // Return 0 if processing fails.
      return y;                         // Return the output sample.
    }                                   // ----------- Read ------------------ //
    // Smooth glide: move to a new delay in k samples.
    void RampTo(
      T newdel,                           // Target delay to ramp to
      size_t k) noexcept                // Number of samples to ramp to the target
    {                                   // ---------- RampTo ----------------- //
      newdel=std::clamp(newdel, T(0.0), T(MaxLen - 1));
      delay=std::clamp(delay, T(0.0), T(MaxLen - 1));
      incr=(newdel-delay)/static_cast<T>(k); // Compute the increment for the ramp.
      targ=newdel;                  // Set the target delay.
      Tick();
    }                                   // ---------- RampTo ----------------- //
    // Call once per sample to progress an active glide.
    void Tick(void) noexcept
    {                                   // ---------- Tick ------------------- //
      // Something happened... update myself.
      if((incr>0&&delay<targ)||(incr<0&&delay>targ))
      {                                 // We updated the incr or delay greater than target?
        delay+=incr;                    // We increased delay by this much.
        // If we are ramping push delay below 1, clamp targ as well
        delay=std::clamp(delay, T(0.0), T(MaxLen - 1));
        mu=delay-std::floor(delay); // Get the fractional part of the delay.
        fd->SetMu(mu); // Update the interpolator with the delay state.
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
      delay=std::clamp(d,T(0.0),T(MaxLen-1)); // Clamp the delay to the range [0, MaxLen-1].
      order=std::min<std::size_t>(N,MaxOrder);// Set the order of the interpolator.
      mu=delay-std::floor(delay); // Get the fractional part of the delay.
      fd->SetOrder(order);               // Set the order of the Lagrange interpolator.
      fdi->SetOrder(order);              // Set the order of the Lagrange deinterpolator.
      fd->SetMu(mu); // Set the fractional part of the delay.
      fdi->SetMu(mu); // Set the fractional part of the delay.
      dl->Clear();                // Clear the delay line buffer.
      incr=T(0.0);                    // Reset the increment to 0.
      targ=delay;            // Set the target delay to the current delay.
      haswritten=false; // Reset the written flag.
    }                                   // ----------- SetDelay ---------------- //
    void SetMu(
      T mu) noexcept                    // Set the fractional part of the delay line
    {                                   // ---------- SetMu ----------------- //
      mu=std::clamp(mu, T(0.0), std::nextafter(T(1.0), T(0.0))); // Clamp the fractional part of the delay line to the range [0, MaxLen-1].
      delay=std::floor(delay)+mu;       // Set the delay to the integer part
      delay=std::clamp(delay,T(0.0),T(MaxLen-1)); // Clamp the delay to the range [0, MaxLen-1].
      targ=delay;
      this->mu=mu;
      fd->SetMu(mu);                    // Set the fractional part of the delay line.
      fdi->SetMu(mu);                   // Set the fractional part of the delay line.
    }          
    void SetIncrement(
      T incr) noexcept                  // Set the increment for the ramp
    {                                   // ---------- SetIncrement ----------------- //
      this->incr=incr;                  // Set the increment for the ramp.
    }                                   // ---------- SetIncrement ----------------- //
    inline void SetOrder(
      size_t order) noexcept             // Set the order of the Lagrange interpolator
    {                                   // ---------- SetOrder ----------------- //
      this->order=std::min<size_t>(order,MaxOrder); // Set the order of the Lagrange interpolator.
      fd->SetOrder(this->order);         // Set the order of the Lagrange interpolator.
    }                                   // ---------- SetOrder ----------------- //
    float PeekIndex(
      size_t idx) const noexcept         // Peek at the sample at the specified index
    {                                   // ---------- PeekIndex ----------------- //
      return dl->PeekIndex(idx);              // Return the sample at the specified index.
    }                                   // ---------- PeekIndex ----------------- //
    inline T GetDelay(void) const noexcept // Get the current delay in samples
    {                                   // ---------- GetDelay ----------------- //
      return delay;                     // Return the current delay in samples.
    }                                   // ---------- GetDelay ----------------- //
    inline T GetMu(void) const noexcept // Get the fractional part of the delay line
    {                                   // ---------- GetMu ----------------- //
      return mu;                       // Return the fractional part of the delay line.
    }                                   // ---------- GetMu ----------------- //
    inline size_t GetOrder(void) const noexcept // Get the order of the Lagrange interpolator
    {                                   // ---------- GetOrder ----------------- //
      return order;                    // Return the order of the Lagrange interpolator.
    }                                   // ---------- GetOrder ----------------- //
    inline size_t GetMaxLen(void) const noexcept // Get the maximum length of the delay line
    {                                   // ---------- GetMaxLen ----------------- //
      return MaxLen;                   // Return the maximum length of the delay line.
    }                                   // ---------- GetMaxLen ----------------- //
    inline size_t GetMaxOrder(void) const noexcept // Get the maximum order of the Lagrange interpolator
    {                                   // ---------- GetMaxOrder ----------------- //
      return MaxOrder;                 // Return the maximum order of the Lagrange interpolator.
    }                                   // ---------- GetMaxOrder ----------------- //
    inline bool HasWritten(void) const noexcept // Check if the delay line has been written to
    {                                   // ---------- HasWritten ----------------- //
      return haswritten;               // Return true if the delay line has been written to, false otherwise.
    }                                   // ---------- HasWritten ----------------- //
    inline void SetSampleRate(
      double fs) noexcept                // Set the sample rate of the delay line
    {                                   // ---------- SetSampleRate ----------------- //
      this->fs=std::max<double>(fs, 1.0); // Set the sample rate of the delay line to the maximum of the provided sample rate and 1.0
      dl->SetSampleRate(fs);            // Set the sample rate of the delay line.
    }                                   // ---------- SetSampleRate ----------------- //
    inline double GetSampleRate(void) const noexcept // Get the sample rate of the delay line
    {                                   // ---------- GetSampleRate ----------------- //
      return fs;                       // Return the sample rate of the delay line.
    }                                   // ---------- GetSampleRate ----------------- //
    private:
      sig::DelayLine<T,MaxLen>* dl{nullptr}; // Delay line buffer
      FarrowInterpolator<T,MaxLen,MaxOrder>* fd{nullptr}; // Farrow filter for
      FarrowDeinterpolator<T,MaxLen,MaxOrder>* fdi{nullptr}; // Deinterpolator for writing samples
      double fs{48000.0}; // Sample rate, default is 44100 Hz.
      size_t order{3};                // Order of the Lagrange filter, default is 3.
      T mu{0.0f};                       // fractional delay in samples
      T delay{T(0)};                     // Current integer delay in samples.
      T incr{0.01};                      // Increment for the ramp.
      T targ{0.0f};  // The target delay                      // Target delay to ramp to.
      bool haswritten{false};
  };                                  // class FarrowDelay
}                                     // namespace sig::wg
