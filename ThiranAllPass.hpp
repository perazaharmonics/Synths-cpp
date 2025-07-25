/*
 * *
 * * Filename: ThiranAllPass.hpp
 * *
 * * Description:
 * *  An N-th-order, maximally-flat, all-pass fractional-delay filter based on
 * *  Thiran’s analytic coefficient set.  Unlike polynomial FIR and Farrow
 * *  structures that trade pass-band ripple for linear phase, the Thiran filter
 * *  achieves **unity magnitude** for every frequency (|H(e^{jw})| = 1) while
 * *  shaping the phase to realise a non-integer group delay of
 * *          D = N - mu , 0 <= mu < 1
 * *  with perfect flatness up to order 2N about ? = 0.  This makes it ideal for:
 * *    • Fine-grain pitch tuning in string and tube waveguides
 * *    • Dispersion simulation without altering amplitude
 * *    • Phase-correct junctions inside feedback delay networks (FDN)
 * *    • Replacing FIR/Lagrange interpolators when modulation depth is small
 * *
 * *  Coefficients:
 * *      a_0 = 1
 * *      a_k = (-1)^k · C(N,k) · PROD_{n=0}^{N-1} (mu - n)/(mu - k - n), 1 <= k <= N
 * *
 * *  Implementation:
 * *   • Direct-Form II Transposed ? O(N) multiplies, O(N) states
 * *   • Header-only template      ? inlinable, no extra linkage
 * *   • Supports float / T / SIMD sample types via the Sample template
 * *
 * * Author:
 * *  JEP  J. Enrique Peraza
 * *
 */
#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <experimental/simd>
#include "DelayLine.hpp"

namespace sig::wg
{
template <typename T=float, size_t MaxLen=1024, size_t MaxOrder=5>  
  class ThiranAllPass
  {
    public:
      ThiranAllPass(void)=default;       // Default constructor
      ~ThiranAllPass(void) noexcept = default; // Default destructor
      bool Prepare(
        size_t order,                   // Order of the Thiran filter, must be >= 1.
        T f) noexcept              // Fractional delay
      {                                 ///%Prepare 
        if (order==0) return false;     // Sanitize input: can't have zero order.
        if (f<0.0||f>=static_cast<T>(order)) 
          return false; // Sanitize input: delay must be in [0, order).
        N=static_cast<int>(order);       // Set the order of the Thiran filter.
        mu=f;                            // Set the fractional delay.
        Assemble(); // Calculate the coefficients for the Thiran filter.
        Clear();                         // Clear the state of the filter.
        return true;                     // Return true if preparation was successful.
      }                                  ///%Prepare 
      inline T ProcessSample(const T& x) noexcept
      {                                 ///%ProcessSample
         auto s=x;                      // Copy the input sample to a temporary variable.
         // --------------------------- //
         // Direct-Form II-Transposed (minimum memory) ladder
         // --------------------------- //
         for (size_t i=0;i<N;++i)       // For each coefficient in Thiran filter
         {                              // Circulate the signal through the filter
            T tmp=s-static_cast<T>(a[i+1])*z[i]; // Compute the intermediate sample using the Thiran coefficients.
            s=z[i]+static_cast<T>(a[i+1])*tmp; // Compute the output sample.
            z[i]=tmp;                    // Update the state register with the new sample.
         }                              // Done circulating the sample
         return s;                      // Return the processed sample.
      }                                 ///%ProcessSample
      inline void Clear(void) noexcept
      {
        std::fill(z.begin(), z.end(), T(0)); // Clear the state registers.
      }
      inline int GetOrder(void) const noexcept
      {                                 //%GetOrder
        return N;                       // Return the order of the Thiran filter.
      }                                 //%GetOrder
      inline T GetFractionalDelay(void) const noexcept
      {                                 //%GetFractionalDelay
        return mu;                     // Return the fractional delay.
      }                                 //%GetFractionalDelay
      inline void SetOrder(int order) noexcept
      {
        if (order < 1 || order>MaxOrder)
          return;
        N=order;
        Assemble();
      }
      inline void SetFractionalDelay(T f) noexcept
      {
        if (f < 0.0 || f >= static_cast<T>(N))
        return;
        mu=f;
        Assemble();
      }
      const std::vector<T>& GetCoefficients(void) const noexcept
      {                                 //%GetCoefficients
        return a;                      // Return the coefficients of the Thiran filter.
      }                                 //%GetCoefficients
    private:
      int N{0};                         // Order of the Thiran filter, must be >= 1.
      T mu{0.0};                   // Fractional delay, 0 <= mu < 1.
      std::vector<T> a;                 // Coefficients of the Thiran filter.
      std::vector<T> z;                 // State registers of the delay line (unit delay).
      // ------- coefficient generation (Thiran 1981, derived from maximally-flat delay) -------
      void Assemble(void) noexcept
      {                                 ///%Assemble:
        a.assign(N+1,0.0);             // Resize the coefficients vector to hold N+1 coefficients.
        a[0]=1.0;                       // First coeff is always 1.0
        // ---------------------------- //
        // For each coeff up to order N, compute the Thiran coefficients
        // ---------------------------- //
        for (int k=1;k<=N;++k)          // 
          a[k]=((k&1)?-1.0:1.0)*Binomial(N,k)*Fracts(k);
        z.assign(N,T{});                // Init the state registers to hold N states.
      }                                 ///%Assemble.
      // ---------------------------- //
      // C(n,k) – symmetrical, numerically stable for small orders (<= 10 normal in audio)
      // ---------------------------- //
      static T Binomial(int n, int k) noexcept
      {                                 //%Binomial
        if (k==0||k==n) return 1.0; // C(n,0) = C(n,n) = 1
        T res=1.0;                 // Initialize the result to 1.0
        for (int i=1;i<=k;++i)          // For each i from 1 to k
          res*=static_cast<T>(n-k+i)/static_cast<T>(i); // Compute the binomial coefficient using the formula C(n,k) = n!/(k!(n-k)!)
        return res;                     // Return the computed binomial coefficient.
      }                                 //%Binomial
      // ---------------------------- //
      // Fracts(k): PROD{n=0}^{N-1} (d - n)/(d - k - n) 
      // see Thiran 1981, eq. 2
      // ---------------------------- //
      T Fracts(int k) const noexcept
      {                                 //%Fracts
        T p=1.0;                   // Initialize the product to 1.0
        for (int n=0;n<N;++n)           // For each n from 0 to N-1
          p*=(mu-n)/(mu-k-n);           // Compute the product of the fractions for the Thiran coefficients.
        return p;                       // Return the computed product.
      }                                 //%Fracts
  };
  template <typename T=float, size_t MaxLen=1024, size_t MaxOrder=5>  
  class ThiranInterpolator
  {

    public:
    bool Prepare(
      float totdel,                     // The total delay in samples.
      size_t order=3) noexcept          // The order of the filter
    {                                   //%Prepare
      // Split integer + fractional
      int D = static_cast<int>(std::floor(totdel));
      float mu = totdel - D;
      this->del = D;                  // Integer delay
      apass.Prepare(order, mu);       // Prepare the Thiran all-pass
      return true;
    }
    inline T Read(size_t i=0) const noexcept
    {
      return apass.GetCoefficients(i);
    }
    // Caller feeds the integer delayed sample here.
    inline T Write(T s) noexcept
    {
      return apass.ProcessSample(s);    // Process sample and return output
    }
    inline T ProcessSample(T x) noexcept
    {
      return apass.ProcessSample(x);     // Process the sample through the Thiran all-pass filter.
    }
    
    inline void SetIntegerDelay(int d) noexcept { del = d; }
    inline void SetFractionalDelay(T f) noexcept
    {
      apass.SetFractionalDelay(f);     // Set the fractional delay in the Thiran all-pass filter.
    }
    inline int GetIntegerDelay() const noexcept { return del; }
    inline T GetFractionalDelay() const noexcept { return apass.GetFractionalDelay(); }
    private:
      int del{0};                       // Integer part of delay in samples
      ThiranAllPass<T> apass;       // Thiran all-pass filter for fractional delay
  }; 
  template <typename T=float,size_t MaxOrder=5>  
  class ThiranDeinterpolator
  {
    public:
      bool Prepare(size_t order, float f) noexcept
      {
        fwd.Prepare(order,f);
        // ---------------------------- //
        // Produce reverse coefficient list
        // ---------------------------- //
        const auto& a=fwd.GetCoefficients(); // Get the coefficients from the forward Thiran filter
        N=static_cast<int>(order);      // Set the order of the Thiran deinterpolator
        b.assign(a.begin(),a.end());    // Copy the coefficients to the reverse list
        z.assign(N,T(0));               // Initialize the state registers to hold N states.
        return true;                    // Return true if preparation was successful.
      }
      inline T ProcessSample(T x) noexcept
      {
        // ---------------------------- //
        // Direct-form II using a (feedback) & b (feed-forward) coefficients
        // ---------------------------- //
        T s=x;                      // Copy the input sample to a temporary variable.
        for (int k=0;k<N;++k)           // For each coefficient in the
        {                               // Circulate samples through filter graph
           T tmp=s-static_cast<T>(ForwardCoeff(k+1)*z[k]);// Difference part of the filter graph
           s=z[k]+static_cast<T>(b[k+1])*tmp;// Compute the output sample using the Thiran coefficients.
           z[k]=tmp;                    // Update the state register with the new sample.
        }                               // Done circulating the sample
        return s;                      // Return the processed sample.
      }                                //%ProcessSample
      inline T Read(size_t i=0) const noexcept
      {
        return fwd.GetCoefficients(i);
      }
      inline void Clear(void) noexcept
      {
        std::fill(z.begin(), z.end(), T(0)); // Clear the state registers.
      }
      inline int GetOrder(void) const noexcept
      {                                 //%GetOrder
        return N;                       // Return the order of the Thiran deinterpolator.
      }                                 //%GetOrder
      inline void SetOrder(int order) noexcept
      {
        if (order < 1 || order>MaxOrder)
          return;                       // Sanitize input: can't have zero order.
        N=order;                        // Set the order of the Thiran deinterpolator.
        b.resize(N+1);                  // Resize the coefficients vector to hold N+1 coefficients
        z.resize(N,T(0));               // Resize the state registers to hold N states.
        fwd.SetOrder(order);            // Set the order of the forward Thiran all-pass filter.
        // ---------------------------- //
        // Produce reverse coefficient list
        // ---------------------------- //
        const auto& a=fwd.GetCoefficients(); // Get the coefficients from the forward Thiran filter
        b.assign(a.begin(),a.end());    // Copy the coefficients to the reverse list
      }                                 //%SetOrder
    private:
      inline T ForwardCoeff(int k) const noexcept { return fwd.GetCoefficients()[k];}
      int N{0};                         // Order of the Thiran deinterpolator
      ThiranAllPass<T> fwd;         // Forward Thiran all-pass filter
      std::vector<T> b;             // Coefficients of the Thiran deinterpolator
      std::vector<T> z;             // State registers of the Thiran de
  };
  template <typename T = float, 
   size_t MaxOrder = 5,
   typename vT = std::experimental::native_simd<T>>
   class ThiranAllPassSIMD
   {
    public:
      static constexpr size_t VL=vT::size();   // Vector length of the SIMD type
      using value_type=T;                      // Value type of the SIMD type
      using packet_type=vT;            // Packet type of the SIMD type
      ThiranAllPassSIMD(void) noexcept
      {
        std::fill(a.begin(), a.end(), T(0)); // Initialize the coefficients to zero.
        std::fill(z.begin(), z.end(), vT(T(0))); // Initialize the state registers to zero.
        stateless=true;                // Set the stateless flag to true.

      }
      ~ThiranAllPassSIMD(void) noexcept = default; // Default destructor
      bool Prepare(                     // Prepare the Thiran all-pass filter
        size_t order,                   // Order of the Thiran filter, must be >= 1.
        T f) noexcept              // Fractional delay
      {                                 //%Process:
        if (order==0||f<0.0||f>=static_cast<T>(order))
          return false;                // Sanitize input: can't have zero order or invalid delay.
        N=static_cast<int>(order);       // Set the order of the Thiran filter
        mu=f;                            // Set the fractional delay.
        Assemble();                      // Calculate the coefficients for the Thiran filter.
        Clear();                         // Clear the state of the filter.
        return true;                     // Return true if preparation was successful.
      }                                  //%Process.
      // Process one  SIMD frame (VL voices)
      inline vT ProcessFrame(vT x) noexcept
      {
        vT s=x;                         // Copy the input sample to a temporary variable.
        // --------------------------- //
        // Direct-Form II-Transposed (minimum memory) ladder
        // --------------------------- //
        for (size_t i=0;i<N;++i)         // For each coefficient in the Thiran filter
        {                                // Circulate the signal through the filter
          vT ak=vT(static_cast<T>(a[i+1])); // Get the Thiran coefficient for the current order
          vT tmp=s-ak*z[i];
          s=z[i]+ak*tmp;                // Compute the output sample using the Thiran coefficients.
          z[i]=tmp;                     // Update the state register with the new sample.
        }
        stateless=false;               // We have circulated through the filter.
        return s;                      // Return the processed sample.
    }                                  //%ProcessFrame
    // Interleaved bufffer helper (framesxVL)
    inline void ProcessBlock(
      const T* src,                     // Inpute buffer
      T* const dst,                     // Output stream
      size_t frames) noexcept           // The frame count
    {                                   //%ProcessBlock
      if (src==nullptr||frames==0) return;
      for (size_t i=0;i<frames;++i)     // For each incoming frame
      {                                 // Break frame in chunks of VL samples
        vT in=vT::load(src+i*VL,std::experimental::parallelism_v2::vector_aligned);// Load the input frame into a SIMD vector
        vT out=ProcessFrame(in);        // Process the frame through the Thiran all
        out.copy_to(dst+i*VL,std::experimental::parallelism_v2::vector_aligned);// Neatly place data block in output buffer.
      }                                // Done processing the block
      stateless=false;                 // We have circulated through the filter.
    }                                  //%ProcessBlock
    inline void Clear(void) noexcept
    {
      for (auto& v: z) v=vT(T(0));
      stateless=true; // Clear the state registers.
    }

    inline int GetOrder(void) const noexcept
    {                                 //%GetOrder
      return N;                       // Return the order of the Thiran filter.
    }                                 //%GetOrder
    inline void SetOrder(int order) noexcept
    {
      if (order < 1 || order>MaxOrder)
        return;                       // Sanitize input: can't have zero order.
      N=order;                        // Set the order of the Thiran filter.
      a.resize(N+1);                  // Resize the coefficients vector to hold N+1 coefficients
      z.resize(N,vT(T(0)));           // Resize the state registers to hold N states.
      Assemble();                     // Assemble the coefficients for the Thiran filter.
      stateless=true;                 // No signal has bothered us yet.
    }                                 //%SetOrder
    inline void SetFractionalDelay(T f) noexcept
    {
      if (f<0.0||f>=static_cast<T>(N))
        return;                       // Sanitize input: can't have zero order.
      mu=f;                            // Set the fractional delay in the Thiran all-pass filter.
      Assemble();                      // Recalculate the coefficients for the Thiran filter.
      stateless=true;                 // No signal has bothered us yet.
    }
    inline T GetFractionalDelay(void) const noexcept
    {                                 //%GetFractionalDelay
      return mu;                     // Return the fractional delay.
    }                                 //%GetFractionalDelay
    inline const std::vector<T>& GetCoefficients(void) const noexcept
    {                                 //%GetCoefficients
      return a;                      // Return the coefficients of the Thiran filter.
    }                                 //%GetCoefficients
    inline bool IsStateless(void) const noexcept
    {
      return stateless;             // True if state registers are cleared, i.e, no input yet.
    }
    private:
      int N{0};                         // Order of the Thiran filter, must be >= 1.
      T mu{0.0};                        // Fractional delay, 0 <= mu < 1.
      std::vector<T> a;                 // Coefficients of the Thiran filter.
      std::vector<vT> z;                 // State registers of the delay line (unit delay).
      bool stateless{true};             // True if unit delays are cleared, i.e, no input yet.
      // ------- coefficient generation (Thiran 1981, derived from maximally-flat delay) -------
      void Assemble(void) noexcept
      {                                 ///%Assemble:
        a.assign(N+1,0.0);             // Resize the coefficients vector to hold N+1 coefficients.
        a[0]=1.0;                       // First coeff is always 1.0
        // ---------------------------- //
        // For each coeff up to order N, compute the Thiran coefficients
        // ---------------------------- //
        for (int k=1;k<=N;++k)          // 
          a[k]=((k&1)?-1.0:1.0)*Binomial(N,k)*Fracts(k);
        z.assign(N,vT(T(0)));                // Init the state registers to hold N states.
      }                                 ///%Assemble.
      // ---------------------------- //
      // C(n,k) – symmetrical, numerically stable for small orders (<= 10 normal in audio)
      // ---------------------------- //
      static T Binomial(int n, int k) noexcept
      {                                 //%Binomial
        if (k==0||k==n) return 1.0; // C(n,0) = C(n,n) = 1
        T res=1.0;                 // Initialize the result to 1.0
        for (int i=1;i<=k;++i)          // For each i from 1 to k
          res*=static_cast<T>(n-k+i)/static_cast<T>(i); // Compute the binomial coefficient using the formula C(n,k) = n!/(k!(n-k)!)
        return res;                     // Return the computed binomial coefficient.
      }                                 //%Binomial
      // ---------------------------- //
      // Fracts(k): PROD{n=0}^{N-1} (d - n)/(d - k - n) 
      // see Thiran 1981, eq. 2
      // ---------------------------- //
      T Fracts(int k) const noexcept
      {                                 //%Fracts
        T p=1.0;                   // Initialize the product to 1.0
        for (int n=0;n<N;++n)           // For each n from 0 to N-1
          p*=(mu-n)/(mu-k-n);           // Compute the product of the fractions for the Thiran coefficients.
        return p;                       // Return the computed product.
      }                                 //%Fracts
   };                                   // %ThiranAllPassSIMD
  template <typename T = float, 
     size_t MaxOrder = 5,
     size_t MaxLen = 1024,
     typename vT = std::experimental::native_simd<T>>
    class ThiranInterpolatorSIMD
  {
    public:
      ThiranInterpolatorSIMD(void) noexcept
      {
        apasses=new ThiranAllPassSIMD<T,MaxOrder,vT>[MaxOrder];  // Our base Thiran all-pass filter array
        lens=new int[MaxOrder];         // Integer delays for each Thiran all-pass filter
        for (size_t i=0;i<MaxOrder;++i) // For the amount of parallel Thiran filters
          lens[i]=0;                    // Initialize the integer delays to zero.
      }
      ~ThiranInterpolatorSIMD(void) noexcept
      {                                 
        delete[] lens;                  // Delete the integer delays array
        lens=nullptr;                   // Remember we deleted it.
        delete[] apasses;               // Delete All-Pass filter bank.
        apasses=nullptr;                // Remember we deleted it.
      }
      bool Prepare( 
        const float* const totdel,                   // The integer + fractional delay in samples.
        const size_t order=3) noexcept        // The order of the filter
      {                                            //%Prepare
        if (totdel==nullptr) return false; // Sanitize input: can't have null pointers.
        if (!apasses) return false;    // Check if the filter bank is initialized
        for (size_t i=0;i<order;++i)   // For each filter in the bank
        {
          // Split integer + fractional delay
          int D=static_cast<int>(std::floor(totdel[i])); // Get the integer part of the delay
          float mu=totdel[i]-D;         // Get the fractional part of the delay
          lens[i]=D; // Set the delay for each filter
        }
        return true;                   // Preparation successful
      }
      inline void WriteFrame(const vT& x) noexcept
      {
        if (!apasses) return;          // Check if the filter bank is initialized
        for (size_t i=0;i<MaxOrder;++i) // For each filter in the bank
          apasses[i].ProcessFrame(x);  // Process the sample through the Thiran all-pass filter.
      }
      inline void WriteBlock(
        const T* src,                     // Input buffer
        T* const dst,                     // Output stream
        size_t frames) noexcept           // The frame count
      {                                   //%WriteBlock
        if (src==nullptr||frames==0) return; // Sanitize input: can't have null pointers or zero frames.
        for (size_t i=0;i<frames;++i)     // For each incoming frame
        {                                 // Break frame in chunks of VL samples
          vT in=vT::load(src+i*vT::size(),std::experimental::parallelism_v2::vector_aligned);// Load the input frame into a SIMD vector
          WriteFrame(in);                // Process the frame through the Thiran all-pass filter.
          in.copy_to(dst+i*vT::size(),std::experimental::parallelism_v2::vector_aligned);// Neatly place data block in output buffer.
        }                                // Done processing the block
      }                                  //%WriteBlock
      inline void SetIntegerDelay(
        size_t i=0, int d=0) noexcept
      {
        if (i>=MaxOrder) return;         // Sanitize input: index must be in range [0, MaxOrder)
        if (d<0) return;                 // Sanitize input: delay must be non-negative
        lens[i]=d;                       // Set the integer delay for the specified filter.
      }
      inline void SetFractionalDelay(
        size_t i=0, float f=0.0f) noexcept
      {
        if (i>=MaxOrder) return;         // Sanitize input: index must be in range [0, MaxOrder)
        if (f<0.0f||f>=static_cast<float>(MaxOrder)) return; // Sanitize input: delay must be in range [0, MaxOrder)
        apasses[i].SetFractionalDelay(f); // Set the fractional delay for the specified filter.
      }
      inline int GetIntegerDelay(size_t i=0) const noexcept
      {
        if (i>=MaxOrder) return 0;       // Sanitize input: index must be in range [0, MaxOrder)
        return lens[i];                  // Return the integer delay for the specified filter.
      }
      inline float GetFractionalDelay(size_t i=0) const noexcept
      {
        if (i>=MaxOrder) return 0.0f;    // Sanitize input: index must be in range [0, MaxOrder)
        return apasses[i].GetFractionalDelay(); // Return the fractional delay for the specified filter.
      }
      inline int GetOrder(void) const noexcept
      {
        return MaxOrder;                // Return the maximum order of the Thiran all-pass filter bank.
      }
      inline const std::vector<T>& GetCoefficients(size_t i=0) const noexcept
      {
        if (i>=MaxOrder) return std::vector<T>{}; // Sanitize input: index must be in range [0, MaxOrder)
        return apasses[i].GetCoefficients(); // Return the coefficients of the specified Thiran all-pass filter.
      }
    private:
      ThiranAllPassSIMD<T,MaxOrder,vT>* apasses{nullptr}; // Forward Thiran all-pass filter
      int* lens{nullptr}; // Integer delays for each Thiran all-pass filter
  };
   template <typename T=float,size_t MaxOrder=5,
    size_t MaxLen = 1024,
    typename vT=std::experimental::native_simd<T>>
    class ThiranDeinterpolatorSIMD
    {
      public:
        ThiranDeinterpolatorSIMD(void) noexcept
        {
          fwd=new ThiranAllPassSIMD<T,MaxOrder,vT>();
          b.resize(MaxOrder, T(0));  // Coefficients of the Thiran deinterpolator
          z.fill(vT(T(0)));  // State
        }
        ~ThiranDeinterpolatorSIMD(void) noexcept
        {
          delete fwd;                // Delete the forward Thiran all-pass filter
          fwd=nullptr;               // Remember we deleted it.
        }
        bool Prepare(
          size_t order,               // Orders of the Thiran deinterpolator, must be >= 1.
          float f) noexcept           // Fractional delay
        {
          fwd->Prepare(order,f);      // Prepare the forward Thiran all-pass filter
          // ------------------------ //
          // Produce reverse coefficient list
          // ------------------------ //
          const auto& a=fwd->GetCoefficients();// Get the coefficients from the forward Thiran filter
          N=static_cast<int>(order);   // Set the order of the Thiran deinterpolator
          b.assign(a.begin(),a.end());  // Copy the coefficients to the reverse list
          z.assign(N,vT(T(0)));         // Initialize the state registers to hold N states.
          return true;                 // Return true if preparation was successful.     
        }
        inline vT ProcessFrame(vT x) noexcept
        {
          // ---------------------------- //
          // Direct-form II using a (feedback) & b (feed-forward) coefficients
          // ---------------------------- //
          vT s=x;                      // Copy the input sample to a temporary variable.
          for (int k=0;k<N;++k)           // For each coefficient in the Thiran deinterpolator
          {                               // Circulate samples through filter graph
            vT ak=vT(ForwardCoeff(k+1));  // Get the Thiran coefficient for the current order
            vT tmp=s-ak*z[k];             // Difference part of the filter graph
            s=z[k]+b[k+1]*tmp;            // Compute the output sample using the Thiran coefficients.
            z[k]=tmp;                     // Update the state register with the new sample.
          }                               // Done circulating the sample
          return s;                      // Return the processed sample.
        }                                 //%ProcessFrame
        // Interleaved bufffer helper (framesxVL)
        inline void ProcessBlock(
          const T* src,                     // Input buffer
          T* const dst,                     // Output stream
          size_t frames) noexcept           // The frame count
        {                                   //%ProcessBlock
          if (src==nullptr||frames==0) return; // Sanitize input: can't have null pointers or zero frames.
          for (size_t i=0;i<frames;++i)     // For each incoming frame
          {                                 // Break frame in chunks of VL samples
            vT in=vT::load(src+i*vT::size(),std::experimental::parallelism_v2::vector_aligned);// Load the input frame into a SIMD vector
            vT out=ProcessFrame(in);        // Process the frame through the Thiran deinterpolator
            out.copy_to(dst+i*vT::size(),std::experimental::parallelism_v2::vector_aligned);// Neatly place data block in output buffer.
          }                                // Done processing the block
        }                                  //%ProcessBlock
        inline void Clear(void) noexcept
        {
          for (auto& v: z) v=vT(T(0)); // Clear the state registers.
        }
        inline int GetOrder(void) const noexcept
        {                                 //%GetOrder
          return N;                       // Return the order of the Thiran deinterpolator.
        }                                 //%GetOrder
        inline void SetOrder(int order) noexcept
        {
          if (order < 1 || order>MaxOrder)
            return;                       // Sanitize input: can't have zero order.
          N=order;                        // Set the order of the Thiran deinterpolator.
          b.resize(N+1);                  // Resize the coefficients vector to hold N+1 coefficients
          z.resize(N,vT(T(0)));           // Resize the state registers to hold N states.
          fwd->SetOrder(order);            // Set the order of the forward Thiran all-pass filter.
          // ---------------------------- //
          // Produce reverse coefficient list
          // ---------------------------- //
          const auto& a=fwd->GetCoefficients(); // Get the coefficients from the forward Thiran filter
          b.assign(a.begin(),a.end());    // Copy the coefficients to the reverse list
        }                                 //%SetOrder
        inline void SetFractionalDelay(float f) noexcept
        {
          if (f<0.0||f>=static_cast<float>(N))
            return;                       // Sanitize input: can't have zero order.
          fwd->SetFractionalDelay(f);      // Set the fractional delay in the Thiran all-pass filter.
          // Recalculate the coefficients for the Thiran deinterpolator.
          const auto& a=fwd->GetCoefficients(); // Get the coefficients from the forward Thiran filter
          b.assign(a.begin(),a.end());    // Copy the coefficients to the reverse list
          z.assign(N,vT(T(0)));           // Initialize the state registers to hold N states.
        }                                 //%SetFractionalDelay
        inline float GetFractionalDelay(void) const noexcept
        {                                 //%GetFractionalDelay
          return fwd->GetFractionalDelay(); // Return the fractional delay.
        }                                 //%GetFractionalDelay
        inline const std::vector<T>& GetCoefficients(void) const noexcept
        {                                 //%GetCoefficients
          return fwd->GetCoefficients(); // Return the coefficients of the Thiran deinterpolator.
        }                                 //%GetCoefficients
      private:
        inline T ForwardCoeff(int k) const noexcept { return fwd->GetCoefficients()[k];}
        int N{0};                // Order of the Thiran deinterpolator
        ThiranAllPassSIMD<T,MaxOrder,vT>* fwd; // Forward Thiran all-pass filter
        std::vector<T> b;      // Thiran inverse coefficients
        std::array<vT,MaxOrder> z; // State registers of the Thiran deinterpolator
      
    };
}
