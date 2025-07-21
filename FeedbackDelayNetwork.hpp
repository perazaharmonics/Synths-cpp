/**
 *  FeedBackDelayNetwork.hpp
 *  -------------------------------------------------------------
 *  A high-quality, compile-time-sized FDN core using FarrowDelayLine
 *  for smooth fractional delays and modulation.  Designed as the
 *  foundation for echo, reverb, chorus/flanger, tremolo, etc.
 *
 *  Template parameters
 *  -------------------
 *    T         : sample type (float / double)
 *    N         : number of delay lines (power of two recommended)
 *
 *  Key Features
 *  ------------
 *  -> Fixed Branch for every line
 *  -> Arbitrary orthogonal / unitary feedback matrix
 *  -> Per-line one-pole damping filters (high-shelf or LPF)
 *  -> Smooth parameter changes (thread-safe)
 *  -> Wet / Dry mix, Output tap-matrix for stereo or multichannel
 *  -> Prepare / Reset pattern identical to JUCE / VST plug-ins
 *  -> Serve as base for echo, reverb, spatial audio effects, and reflection simulations.
 *
 *  Dependencies: FarrowDelayLine.hpp, FilterFactory.hpp
 */
#pragma once
#include <vector>
#include <array>
#include <atomic>
#include <cmath>
#include <random>
#include <cstddef>
#include <memory>
#include "FilterFactory.hpp"
#include "OnePole.hpp"
#include "BiQuad.hpp"
#include "DelayBranch.hpp"   // use DelayBranch for fractional delays
#include "spectral/Matrices.hpp"
#include "spectral/MatrixSolver.hpp"

namespace sig::wg {
   
  namespace detail                        // Helper namespace 
  {
    template <typename T, std::size_t N> // N filter taps.
    constexpr std::array<std::array<T, N>, N> identity() noexcept
    {                                   // Build Identity matrix (itself?)
        std::array<std::array<T, N>, N> mat{}; // Create an N x N matrix
        for (size_t i=0;i<N;++i)        // Loop through each row
          for (size_t j=0;j<N;++j)      // Loop through each column
            mat[i][j] = (i == j)?T(1):T(0); // Set diagonal elements to 1, others to 0
        return mat;                     // Return the identity matrix
    }

    // In-place householder matrix for maximum echo density
    template <typename T, size_t N>
    constexpr std::array<std::array<T, N>, N> householder(void) noexcept
    {
      static_assert(N>=2,"Householder matrix requires N >= 2");
      constexpr T s=static_cast<T>(2)/static_cast<T>(N);
      std::array<std::array<T, N>,N> m{};
      for (size_t i=0;i<N;++i)
        for (size_t j=0;j<N;++j)
          m[i][j]=(i==j)?s-T(1):s;
      return m;
    }

    template <typename T, std::size_t N>
    constexpr auto hadamard() noexcept
    {
        static_assert((N & (N - 1)) == 0, "Hadamard requires power-of-two N");
        std::array<std::array<T, N>, N> m{};
        m[0][0] = static_cast<T>(1);
        for (std::size_t size = 1; size < N; size <<= 1)
          for (std::size_t i = 0; i < size; ++i)
            for (std::size_t j = 0; j < size; ++j)
            {
              m[i+size][j]= m[i][j];
              m[i][j+size]=m[i][j];
              m[i+size][j+size]=-m[i][j];
            }
        const T norm = static_cast<T>(1)/std::sqrt(static_cast<T>(N));
        for (auto& row:m) for (auto& v:row) v*=norm;
        return m;
    }

    // Random orthogonal (QR of random matrix, seeded once)
    template <typename T, std::size_t N>
    inline auto randomOrthogonal(unsigned seed = 0xCEFAFD) noexcept
    {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<T> dist;
        std::array<std::array<T, N>, N> a{};
        for (auto& row : a)
            for (auto& v : row) v = dist(gen);

        // Gram–Schmidt QR → Q
        for (std::size_t k = 0; k < N; ++k)
        {
          for (std::size_t i = 0; i < k; ++i)
          {
            T dot{};
            for (std::size_t n = 0; n < N; ++n) dot += a[k][n] * a[i][n];
            for (std::size_t n = 0; n < N; ++n) a[k][n] -= dot * a[i][n];
          }
          T norm{};
          for (auto v : a[k]) norm += v * v;
          norm = std::sqrt(norm);
          for (auto& v : a[k]) v /= norm;
        }
        return a;
    }

    enum MatrixKind { Identity, Householder, Hadamard, RandomOrthogonal };

    template <typename T, size_t N>
    inline auto MakeMatrix(MatrixKind k, unsigned seed = 0) noexcept
    {
        switch (k)
        {
        case MatrixKind::Identity:      return identity<T, N>();
        case MatrixKind::Householder:   return householder<T, N>();
        case MatrixKind::Hadamard:      return hadamard<T, N>();
        default:                        return randomOrthogonal<T, N>(seed);
        }
    }

    template<typename T, std::size_t N>
    constexpr std::array<std::array<T,N>,N> toStdArray(const Matrices<T>& M)
    {
        if (M.Rows() != N || M.Cols() != N)
            throw std::invalid_argument{"Matrix size ≠ FDN order"};

        std::array<std::array<T,N>,N> out{};
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                out[i][j] = M(i,j); // assumes Matrices<T>::operator()(r,c)
        return out;
    }     
  } // namespace detail

  // =====================================================================
  // Feedback Delay Network
  // =====================================================================  
  template <typename T=float,            // Our sample format
            size_t MaxLen=1<<15,        // Maximum length of the delay line in samples
            size_t K=3,                 // Farrow filter order
            size_t P=3,                 // Thiran order 
            size_t Ntaps=8>             // Number of delay lines (taps)
  class FeedbackDelayNetwork
  {
      using Branch = DelayBranch<T,MaxLen,K,P>;
  public:
      FeedbackDelayNetwork(void) noexcept
      {
          fbmtx = sig::wg::detail::householder<T, Ntaps>(); // Init feedback matrix
          mtxk  = sig::wg::detail::MatrixKind::Householder; // Track kind
          for (size_t i=0;i<Ntaps;++i)
              dls[i].Prepare(static_cast<size_t>(defaultDelays[i] * this->fs), 0.0f); // Use Prepare instead of SetDelay
      }
      ~FeedbackDelayNetwork(void) noexcept = default;

      // Prepare the FDN with a given delay time and damping factor
      bool Prepare(double fs,size_t bs)
      {                                 // Prepare the FDN with a given sample rate and block size
        if (fs<=0.0||bs<1) return false; // If sample
        this->fs=fs;                    // Set the sample rate
        this->bs=bs;                    // Set the block size
        if (predel>0.0)                 // If we have a predelay
          predelay.Prepare(static_cast<size_t>(predel*fs),0.0f);
        for (size_t i=0;i<Ntaps;++i)    // For each tap (matrix row)
          shelf[i]=filterFactory.LowShelf(fs,shfc,shboost,slope);
        for (size_t i=0;i<Ntaps;++i)    // For each tap, 
        {                               // Set damper filter
          dampLP[i]=filterFactory.OnePoleLP(fs,dfc);// This filter from the factory.
          dampLP[i].Prepare(fs,bs);     // Prepare the damper filter
        }                               // Done preparing the FDN.
        Clear();                        // Clear the FDN state
        return true;                    // Return true if preparation was successful
      }                                 // Prepare the FDN with a given delay time and damping factor

      // Process a block of samples through the FDN
      bool Process(
        const T* in,                    // Input buffer (nullptr for dry output)
        T* const outL,                  // Left output buffer pointer to stream
        T* const outR,                  // Right output buffer pointer to stream
        size_t M) noexcept              // Process a block of samples through the FDN 
      {                                 // ----------- Process --------------------
        for (size_t n=0; n<M;++n)    // For each sample in the block
        {
        bool nofb=true;             // Record if feedback state is empty
        for (size_t i=0;i<Ntaps&&nofb;++i)// For each tap (line)
          for (size_t j=0;j<Ntaps;++j) // For each tap (column)
            if (fbmtx[i][j]!=T(0)) nofb=false; // For the feedback matrix....
        T x=T(0);                   // Initialize output sample
        if (in!=nullptr)            // Is the input null?
        {                           // No, we can go on
          if (wetMix.load(std::memory_order_relaxed)>= T(1.0)&&nofb)
            x=in[n];                // Copy the input sample
          if (predel<=0.0)          // Did we set a predelay?
            x=in[n];                // Yes just copy the input buffer.
          else                      // Else we did set a predelay
          {                         // So configure it.
            predelay.Write(in[n]);  // Write into predelay
            predelay.Propagate(1);  // Circulate ~~~~~~~
            x=predelay.Read();      // Read the predelay output
          }
        }
        const T wet=wetMix.load(std::memory_order_relaxed);
        const T dry=static_cast<T>(1)-wet;
        // If no feedback and wet=1.0, just copy input to output
        if (in && wet==T(1.0)&&nofb&&predel<=T(0.0))
        {                             // ~~~~~~~~~~~~~~~~~~~~~~
          outL[n]=x;                  // Laaaaaaazy
          outR[n]=x;                  //          Coooooooopy
          continue;                   // Outpuuuuut
        }                             // Done doing lazy copy.
        // Gather feedback outputs
        for (size_t i=0;i<Ntaps;++i)  // For the number of coeffs in filer
        {                             // Feed into the filter blocks
          T r = dls[i].Read();      // Read the delay line output
          // ----------------------- //
          // Apply shelf and damping filters
          // ----------------------- //
          T d=(dampfc[i]<0.5*fs)?dampLP[i].ProcessSample(r) : r;
          // Apply the shelf and damper only if below Nyquist limit!
          lastOut[i] = (shelffc[i]<0.5*fs)?shelf[i].ProcessSample(d):d;
        }                             // Done applying filtersssz.
        // -------------------------- //
        // Mix through feedback matrix
        // -------------------------- //
        std::array<T,Ntaps> feed{};   // Initialize feedback array
        for (size_t i=0;i<Ntaps;++i)  // For each tap (rows)
          for (size_t j=0;j<Ntaps;++j)// For each tap (columns)
            feed[i]+=fbmtx[i][j]*lastOut[j];// Mixed output
        // -------------------------- //
        // Data dumping into each delay line
        // -------------------------- //
        for (size_t i=0;i<Ntaps;++i)  // For each tap (row)
        {
          dls[i].Write(x+feed[i]);    // Write the input + feedback to the delay line
          dls[i].Propagate(1);        // Make it Tick().....
        }
          // -------------------------- //
        // Simple stereo tap: even ➜ L, odd ➜ R
        // -------------------------- //
        T yL{},yR{};                  // Initialize left and right outputs
        for (size_t i=0;i<Ntaps;++i) ((i&1)?yR:yL) += lastOut[i];
        const T norm = static_cast<T>(2)/static_cast<T>(Ntaps);
        yL*=norm;yR*=norm;            // Normalize the outputs
        if (in!=nullptr)              // Do we even have an input?
        {                             // Yeah ok let's mix in.
          outL[n]=dry*in[n]+wet*yL;   // Calculate left output
          outR[n]=dry*in[n]+wet*yR;   // Calculate right output
        }                             // Else we don't have an input
        else                          // So just output the wet signal
        {                             // ~~~~~~~~~~~~~~~~~~~~~~~~~~      
          outL[n]=wet*yL;             // Output left wet signal
          outR[n]=wet*yR;             // Output right wet signal
        }                             // Done mixing the output.
       }                               // End of processing block
       return true;                    // Return true if processing was successful
      }                                 // Process a block of samples through the FDN
      void Clear(void) noexcept         // Reset the FDN state
      {                                 // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (auto& sh:shelf) sh.Reset();// Cleeaarr
        predelay.Clear();             //             aaalll
        for (auto& dl:dls) dl.Clear();//                 ddddeelllaaaaaaysssss
        for (auto& d:dampLP) d.Reset();//         and filters.....
        std::fill(lastOut.begin(),lastOut.end(),T(0));// Reset last output
      }                                 // ----------- Clear ------------- //
      // Wet mix utilities
      void SetWetMix(T mix) noexcept
      {
          if (mix<static_cast<T>(0))      mix = static_cast<T>(0);
          else if (mix>static_cast<T>(1)) mix = static_cast<T>(1);
          wetMix.store(mix,std::memory_order_relaxed);
      }
      T GetWetMix(void) const noexcept { return wetMix.load(std::memory_order_relaxed); }

      // Feedback matrix utilities
      void SetFeedbackMatrix(detail::MatrixKind k, unsigned seed=0) noexcept
      {
          fbmtx = detail::MakeMatrix<T,Ntaps>(k,seed);
          mtxk  = k;
      }
      void SetFeedbackMatrix(const std::array<std::array<T,Ntaps>,Ntaps>& m) noexcept
      {
          fbmtx = m;
          mtxk  = detail::MatrixKind::Identity;
      }
      const std::array<std::array<T,Ntaps>,Ntaps>& GetFeedbackMatrix() const noexcept { return fbmtx; }

      using Matrices = std::vector<std::vector<T>>;
      void SetFeedBackMatrix(const Matrices m) { fbmtx = detail::toStdArray<T,Ntaps>(m); }

      // Delay helpers
      void SetPreDelaySeconds(double secs) noexcept
      {
          predel = secs;
          predelay.Prepare(static_cast<size_t>(predel * fs), 0.0f); // Use Prepare instead of SetDelay
      }
      // --------------------------------------------------------------------- //
      //  Compatibility shims – keep old client code alive with the new engine //
      // --------------------------------------------------------------------- //
      inline void SetDelay(double samples) noexcept
      {
        const auto N=static_cast<size_t>(std::floor(samples));
        const auto mu=static_cast<float>(samples-static_cast<double>(N));
        for (auto& dl :dls)
          dl.Prepare(N, mu); // Use Prepare to set delay and fractional delay
      }

      inline void Tick(void) noexcept
      {
        for (auto& dl:dls)      // For # of delay branches in system
          dl.Propagate(1);        // Make it Tick().....
      }
      const std::array<Branch,Ntaps>& GetDelayLines(void) const noexcept { return dls; }
      void SetDelaySeconds(const std::array<double,Ntaps>& secs) noexcept
      { for (size_t i=0;i<Ntaps;++i) dls[i].SetDelay(secs[i]*this->fs); }
      void SetFractionalDelay(const std::array<double,Ntaps>& secs) noexcept
      { for (size_t i=0;i<Ntaps;++i) dls[i].SetFractionalDelay(secs[i]*this->fs); }

      // Damper & shelf helpers
      void SetDamperCutoffs(const std::array<double,Ntaps>& freq) noexcept
      {
        dampfc=freq;                     // Set the cutoff frequencies
        for (size_t i=0;i<Ntaps;++i)     // For each tap (line)
        {                                //     of the FDN
          dampLP[i].Reset();             // Reset the filter
          if (freq[i]<fs*0.5)            // If below Nyquist
          {                              // activate the filter
            dampLP[i] = filterFactory.OnePoleLP(fs,freq[i]);
            dampLP[i].Prepare(fs,bs);    // and assemble it.
          }                              // Done setting the damper cutoffs
        }                                // 
      }
      void SetDamperCutoff(size_t idx, double freq) noexcept
      {
          if (idx>=Ntaps) return;
          dampfc[idx]=freq;
          if (freq<fs*0.5)
              dampLP[idx] = filterFactory.OnePoleLP(fs,freq);
          dampLP[idx].Prepare(fs,freq);
      }

      void SetShelfCutoffs(const std::array<double,Ntaps>& freq) noexcept
      {
        shelffc=freq;
        for (size_t i=0;i<Ntaps;++i)
          if (freq[i]<fs*0.5) shelf[i].SetCutoffFrequency(freq[i]);
      }
      void SetShelfCutoff(int idx, double freq) noexcept
      {
          if (idx<0 || idx>=static_cast<int>(Ntaps)) return;
          shelffc[idx]=freq;
          if (freq<fs*0.5) shelf[idx].SetCutoffFrequency(freq);
      }
      void SetShelfSlopes(const std::array<double,Ntaps>& slope) noexcept
      { for (size_t i=0;i<Ntaps;++i) shelf[i] = filterFactory.LowShelf(fs,shfc,shboost,slope); }
      void SetShelfSlope(int idx, double slope) noexcept
      { if (idx>=0 && idx<static_cast<int>(Ntaps)) shelf[idx] = filterFactory.LowShelf(fs,shfc,shboost,slope); }

      // Delay getters
      inline std::array<double,Ntaps> GetDelays(void) const noexcept
      {
          std::array<double,Ntaps> delays{};
          for (size_t i=0;i<Ntaps;++i) delays[i]=dls[i].GetDelay();
          return delays;
      }
      inline void SetDelays(const std::array<double,Ntaps>& delays) noexcept
      { for (size_t i=0;i<Ntaps;++i) dls[i].Prepare(static_cast<size_t>(delays[i]*this->fs), 0.0f); }

  private:
      double fs{48000.0};
      size_t bs{256};
      Branch predelay;
      double predel{0.3};        // Pre-delay time in seconds
      std::array<Branch,Ntaps> dls;
      std::array<double,Ntaps> defaultDelays{ 0.0,0.0,0.0 };
      std::array<BiQuad<T>,Ntaps> shelf;
      FilterFactory<float> filterFactory;
      std::array<OnePole<T>,Ntaps> dampLP;
      std::array<double,Ntaps>  dampfc{};
      std::array<double,Ntaps> shelffc{};
      double dfc{fs*0.5};
      double shfc{500.0};
      double shboost{0.0};
      double slope{1.0};
      std::array<std::array<T,Ntaps>,Ntaps> fbmtx;
      std::array<T,Ntaps>     lastOut{};
      std::atomic<T>          wetMix{static_cast<T>(0.5)};
      detail::MatrixKind      mtxk = detail::Identity;

      static void Normalize(std::array<std::array<T,Ntaps>,Ntaps>& M) noexcept
      {
          T maxv=T(0);
          for (const auto& row:M) for (const auto& v:row) maxv = std::max(maxv,std::fabs(v));
          if (maxv!=T(0)) for (auto& row:M) for (auto& v:row) v/=maxv;
      }
  };
} // namespace sig::wg
