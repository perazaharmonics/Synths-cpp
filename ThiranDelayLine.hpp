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
      int D=static_cast<int>(std::floor(totdel));
      float mu=totdel-D;
      this->del=D;                    // caller must implement / own int delay!
      apass.Prepare(order,mu);        // Prepare the Thiran all-pass
      return true;                    // Return true if preparation was successful.
    }                                 //%Prepare
    
    // Caller feeds the integer delayed sample here.
    inline void Write(float s) noexcept
    {
      apass.ProcessSample(s);          // Process the sample through the Thiran all-pass filter.
    }
    inline void SetIntegerDelay(int d) noexcept
    {
      del=d;                           // Set the integer part of the delay in samples.
    }
    inline void SetFractionalDelay(float f) noexcept
    {
      apass.SetFractionalDelay(f);     // Set the fractional delay in the Thiran all-pass filter.
    }
    int GetIntegerDelay(void) const noexcept
    {
      return del;                      // Return the integer part of the delay in samples.
    }
    T GetFractionalDelay(void) const noexcept
    {
      return apass.GetFractionalDelay();
    }
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
   template <typename T=float,
     size_t MaxLen=1024,
     size_t MaxOrder=5>
  class ThiranDelayLine
  {
    public:
      ThiranDelayLine(void) noexcept
      {
        dl=new sig::DelayLine<T,MaxLen>(); // Create a new delay line with the specified maximum length.
        apass=new ThiranAllPass<T>();
        inv=new ThiranDeinterpolator<T>();
        Prepare(T(0));
      }
      ~ThiranDelayLine(void)
      {
        delete apass;
        apass=nullptr;
        delete dl;
        dl=nullptr;
        delete inv;
        inv=nullptr;
      }
      void Prepare(
        T d,                         // Total delay int and frac part
        size_t o=3) noexcept         // Order of Thiran All-Pass filters
      {
        order=o;
        delay=d;
        SetDelay(delay,order);
        dl->Clear();                  // Clear the delay line.
      }
      // Hard-set exact delay (no ramp)
      void SetDelay(T d, size_t o=3) noexcept
      {
        delay=std::clamp(d,T(0),T(MaxLen-1)); // Clamp the delay to the range [0, MaxLen-1]
        order=std::clamp<std::size_t>(o,1,MaxOrder);
        Assemble();                 // Assemble the Thiran all-pass filter with the specified order and fractional delay.
      }
      // Smooth glide
      void RampTo(T d,size_t k=1) noexcept
      {
        targ=std::clamp(d,T(0),T(MaxLen-1));
        incr=(targ-delay)/T(std::max<std::size_t>(k,1));
      }
      inline void Tick(void) noexcept
      {
        if ((incr>0&&delay<targ)||(incr<0&&delay>targ))
        {
          delay+=incr;               // We updated the delay by this much.
          delay=std::clamp(delay,T(0),T(MaxLen-1)); // Clamp the delay to the range [0, MaxLen-1]
          Assemble();                // Reassemble the Thiran all-pass filter with the new delay.
        }
      }
      inline void Write(T x) noexcept
      {
        dl->Write(inv.ProcessSample(x));
      }
      inline T Read(void) noexcept
      {
        return ApplyReadFilter(idelay,mu); // Read the sample from the delay line and apply the Thiran all-pass filter.
      }
      // Arbitrary fractional access (forward)
      inline T ReadFrac(T d) noexcept
      {
        int D=int(std::floor(d));
        T mu=T(d-T(D));
        return ApplyReadFilter(D,mu);
      }
      // Tap *backwards* (useful for reverse propagation branches)
      inline T ReadReverse(T d) noexcept
      {
        int D=idelay-int(std::floor(d));
        if (D<0) return T(0); // Sanitize input: D must be in the range [0, MaxLen)
        T mu=this->mu+T(d-std::floor(d));
        if (mu>=1.0)
        {
          ++D;
          mu-=1.0;
        }
        return ApplyReadFilter(D,mu); // Read the sample from the delay line and apply the Thiran all-pass filter.
      }
      // Accessor's and mutators
      inline T GetDelay(void) const noexcept
      {
        return delay;              // Return the current delay in samples.
      }
      inline int GetIntegerDelay(void) const noexcept
      {
        return idelay;            // Return the integer part of the delay in samples.
      }
      inline T GetFractionalDelay(void) const noexcept
      {
        return mu;                // Return the fractional part of the delay.
      }
      inline size_t GetOrder(void) const noexcept
      {
        return order;             // Return the order of the Thiran all-pass filter.
      }
      inline void SetOrder(size_t o) noexcept
      {
        if (o < 1 || o > MaxOrder)
          return;                 // Sanitize input: can't have zero order.
        order=o;                  // Set the order of the Thiran all-pass filter.
        Assemble();               // Reassemble the Thiran all-pass filter with the new order.
      }
      inline void SetFractionalDelay(T f) noexcept
      {
        mu=std::clamp(f,T(0),T(1)); // Clamp the fractional delay to the range [0, 1]
        Assemble();                 // Reassemble the Thiran all-pass filter with the new fractional delay.
      }
    private:
      inline void Assemble(void) noexcept
      {
        idelay=int(std::floor(delay));
        mu=T(delay-T(idelay));
        apass->Prepare(order,mu);
        inv->Prepare(order,mu);
      }
      inline T ApplyReadFilter(int D,T m) noexcept
      {
        if (D<0||D>=int(MaxLen)) return T(0);
        if (std::abs(m)<std::numeric_limits<T>::epsilon())
          return dl->Peek(D); // No fractional delay, just return the sample at D
        ThiranAllPass<T>* tmp=new ThiranAllPass<T>{}; // Create a new Thiran all-pass filter
        tmp->Prepare(order,m); // Prepare the Thiran all-pass filter with the fractional delay
        T x=tmp->ProcessSample(dl->Peek(D)); // Process the sample through the Thiran
        delete tmp; // Clean up the temporary Thiran all-pass filter
        tmp=nullptr; // Set the pointer to null
        return x; // Return the processed sample
      }
    private:
      sig::DelayLine<T,MaxLen>* dl{nullptr};
      ThiranAllPass<T>* apass{nullptr};
      ThiranDeinterpolator<T>* inv{nullptr};
      // state
      T delay{0};                 // Delay with integer and frac part.
      int idelay{0};              // Integer part of the delay in samples
      T incr{0};                 // Increment for the delay
      T targ{0};                 // Target delay in samples
      T mu{0};                   // Fractional part of the delay, default is 0
      size_t order{3};          // Order of the Thiran filter, default is 3
  };                           // Thiran Delay Line
  template <typename T=float,
    size_t MaxLen=1024,
    size_t MaxOrder=5,
    typename packet=std::experimental::native_simd<T>>
  class ThiranDelayLineSIMD
  {
    static constexpr size_t VL=packet::size(); // Vector length of the SIMD type
    public:
      ThiranDelayLineSIMD(void) noexcept
      {
        dl=new sig::DelayLineSIMD<T,MaxLen,packet>{};
        for (size_t i=0;i<VL;++i)       // For the number of voices (Signals)
        {
          inter[i]=new ThiranAllPassSIMD<T,MaxOrder,packet>{}; // Create a new Thiran all-pass filter for each voice.
          deinter[i]=new ThiranAllPassSIMD<T,MaxOrder,packet>{}; // Create a new Thiran all-pass filter for each voice.
          inter[i]->SetOrder(order);     // Set the order of the Thiran all-pass filter.
          deinter[i]->SetOrder(order);   // Set the order of the Thiran all-pass filter.
          inter[i]->SetMu(mu);           // Set the fractional part of the delay line.
          deinter[i]->SetMu(mu);         // Set the fractional part of the delay line.
        }
      }
      ~ThiranDelayLineSIMD(void)
      {
        delete dl;                      // Delete the delay line.
        dl=nullptr;                     // Set the delay line pointer to null.
        for (size_t i=0;i<VL;++i)       // For the number of voices (Signals)
        {
          delete inter[i];              // Delete the Farrow interpolator.
          inter[i]=nullptr;             // Set the Farrow interpolator pointer to null.
          delete deinter[i];            // Delete the Farrow deinterpolator.
          deinter[i]=nullptr;           // Set the Farrow deinterpolator pointer to null.
        }
      }
      void Prepare(T d, size_t o=3) noexcept
      {
        for (int i=0;i<VL;++i) // For each voice (signal)
        {
          inter[i]->Prepare(o,d);       // Prepare the Thiran all-pass filter with the specified order and fractional delay.
          deinter[i]->Prepare(o,d);     // Prepare the Thiran all-pass filter with the specified order and fractional delay.
        }
      }
      void SetDelay(
        T d,                             // The integer and frac delay
        size_t o=3) noexcept             // The order of the MF FIR filters
      {
        delay=std::clamp(d,T(0),T(MaxLen-1));
        order=std::min<size_t>(o,MaxOrder); // Set the order of the Lagrange interpolator.
        Assemble();                      // Assemble the MF FIR filters.
      }
      inline void SetMu(T m) noexcept
      {
        mu=std::clamp(m,T(0),std::nextafter(T(1),T(0))); // Clamp the fractional part of the delay line to the range [0, MaxLen-1].
        delay=std::floor(delay)+mu;     // Set the delay to the integer part plus the fractional part.
        delay=std::clamp(delay,T(0),T(MaxLen-1)); // Clamp the delay to the range [0, MaxLen-1].
        Assemble();                      // Reassemble the MF FIR filters with the new fractional delay.
      }
      void RampTo(
        T d,                             // The int and frac delay
        size_t k) noexcept               // Number of samples to ramp to target
      {
        targ=std::clamp(d,T(0),T(MaxLen-1)); // Clamp the target delay to the range [0, MaxLen-1].
        incr=(targ-delay)/T(std::max<size_t>(k,1));
      }
      inline void Tick(void) noexcept
      {
        if ((incr>0&&delay<targ)||(incr<0&&delay>targ))
        {                                 // We updated the incr or delay greater than target?
          delay+=incr;                    // We increased delay by this much.
          delay=std::clamp(delay,T(0),T(MaxLen-1)); // Clamp the delay to the range [0, MaxLen-1].
          Assemble();                     // Assemble the MF FIR filters with the new delay.
        }                                 // Done updating out state.
      }
      inline void Write(packet x) noexcept
      {
        delay=std::clamp(delay,T(0),T(MaxLen-1)); // Clamp the delay to the range [0, MaxLen-1].
         for (size_t i=0;i<VL;++i)
            x[i]=deinter[i]->ProcessSample(x[i]);
        dl->WriteAt(0,x);
        dl->Advance();                  // circularly move write head.
        haswritten=true; // Mark that we have written a sample.
      }
      //      const DelayLine<T,MaxLen>& dl,   // The delay line to process
      //size_t D,                     // The fractional delay in samples
      //T* const y
      inline packet Read(void) noexcept
      {
        packet y{};                    // Output sample
        for (size_t i=0;i<VL;++i) // For each voice (signal)
        {
          inter[i]->Prepare(order,mu); // Prepare the Thiran all-pass filter with
          T yi=inter[i]->ProcessSample(dl->PeekScalar(idelay,i));
          y[i]=yi;                     // Store the output sample in the output vector.
        }
        return y;                      // Return the output vector containing the output samples for each voice.
      }
      packet ReadFrac(T d) noexcept
      {
        int D=int(std::floor(d)); // Get the integer part of the delay.
        T m=d-T(D);               // Get the fractional part of the delay.
        packet y{};                 // Output sample
        for (size_t i=0;i<VL;++i) // For each voice (signal)
        {
          inter[i]->Prepare(order,m); // Prepare the Thiran all-pass filter with
          T yi=inter[i]->ProcessSample(dl->PeekScalar(idelay,i));
          y[i]=yi;                     // Store the output sample in the output vector.
        }
        return y;                  // Return the output vector containing the output samples for each voice.
      }
    private:
      sig::DelayLineSIMD<T,MaxLen,packet>* dl{nullptr};
      std::array<ThiranAllPassSIMD<T,MaxOrder,packet>,VL> inter{};
      std::array<ThiranAllPassSIMD<T,MaxOrder,packet>,VL> deinter{};
      T fs{48000.0}; // Sample rate, default is 44100 Hz
      size_t order{3};                // Order of the Lagrange filter, default is 3.
      T mu{0.0f};                       // fractional delay in samples
      T delay{T(0)};                     // Current integer delay in samples.
      int idelay{0}; // Current integer delay in samples.
      T incr{0.01};                      // Increment for the ramp.
      T targ{0.0f};  // The target delay                      // Target delay to ramp to.
      bool haswritten{false}; // Flag to indicate if the delay line has been
      inline void Assemble(void) noexcept
      {
        idelay=int(std::floor(delay));
        mu=delay-T(idelay); // Get the fractional part of the delay.
        for (size_t i=0;i<VL;++i) // For each voice (signal)
        {
          inter[i]->Prepare(order,mu); // Prepare the Thiran all-pass filter with the specified order and fractional delay.
          deinter[i]->Prepare(order,mu); // Prepare the Thiran deinterpolator
        }                               // Done with the coefficients.
      }
  };  
}
