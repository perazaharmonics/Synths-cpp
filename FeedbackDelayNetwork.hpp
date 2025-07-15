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
 *  -> Fixed FarrowDelayLine<T> for every line
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
#include "FarrowDelayLine.hpp"
#include "FilterFactory.hpp"
#include "OnePole.hpp"
#include "BiQuad.hpp"
#include "spectral/Matrices.hpp"
#include "spectral/MatrixSolver.hpp"
namespace sig::wg
{
  template <typename T>
  namespace detail
  {
    // Build a NxN identity matrix
    template <typename T, size_t N>
    constexpr std::array<std::array<T,N>,N> identity() noexcept
    {
        std::array<std::array<T,N>,N> mat{}; // Create a NxN matrix
        for (size_t i=0;i<N;++i)             // Loop through each row
            for (size_t j=0;j<N;++j)           // Loop through each column
            mat[i][j]=(i==j)?T(1):T(0);       // Set diagonal elements to 1, others to 0
        return mat;                          // Return the identity matrix
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
                    m[i + size][j]       =  m[i][j];
                    m[i][j + size]       =  m[i][j];
                    m[i + size][j + size]= -m[i][j];
                }
        const T norm = static_cast<T>(1) / std::sqrt(static_cast<T>(N));
        for (auto& row : m) for (auto& v : row) v *= norm;
        return m;
    }
    /* Random orthogonal (QR of random matrix, seeded once)               */
    template <typename T, std::size_t N>
    inline auto randomOrthogonal(unsigned seed = 0xCEFAFD) noexcept
    {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<T> dist;
        std::array<std::array<T, N>, N> a{};
        for (auto& row : a)
            for (auto& v : row) v = dist(gen);

        // Gram–Schmidt QR ? Q
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
    enum MatrixKind
    {
        Identity,
        Householder,
        Hadamard,
        RandomOrthogonal
    };
    template <typename T, size_t N>
    inline auto MakeMatrix(MatrixKind k, unsigned seed = 0) noexcept
    {
        switch (k)
        {
        case MatrixKind::Identity:          return identity<T, N>();
        case MatrixKind::Householder:       return householder<T, N>();
        case MatrixKind::Hadamard:          return hadamard<T, N>();
        default:                            return randomOrthogonal<T, N>(seed);
        }
    }
   template<typename T, std::size_t N>
    constexpr std::array<std::array<T,N>,N>
    toStdArray(const Matrices<T>& M)
    {
        if (M.Rows() != N || M.Cols() != N)
            throw std::invalid_argument{"Matrix size ? FDN order"};

        std::array<std::array<T,N>,N> out{};
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                out[i][j] = M(i,j);         // assumes Matrices<T>::operator()(r,c)
        return out;
    }     
  } // namespace detail
    template <typename T=float, size_t Ntaps=8>
    class FeedbackDelayNetwork
    {
      public:
        FeedbackDelayNetwork(void) noexcept
        {
          fbmtx=detail::householder<T,Ntaps>(); // Initialize the feedback matrix to a householder matrix
        }
        ~FeedbackDelayNetwork(void) noexcept = default; // Default destructor
        // Prepare the FDN with a given delay time and damping factor
        bool Prepare(double fs,size_t bs)
        {
          fs=fs;
          bs=bs;
          predelay.SetDelay(predel); // Set the pre-delay line with the given delay
          for (size_t i=0;i<Ntaps;++i)
            shelf[i]=filterFactory.LowShelf(fs,shfc,shboost,slope); // Create low-shelf filters for each delay line
          for (size_t i=0;i<Ntaps;++i)
            dls[i].SetDelay(predel);     // Prepare the delay line with the given sample rate
          for (size_t i=0;i<Ntaps;++i)
          {
            dampLP[i]=filterFactory.OnePoleLP(fs,dfc); // Create one-pole low-pass filters for damping
            dampLP[i].Prepare(fs,bs);       // Prepare the damping filter with the given sample rate
          }
          Clear();        // Clear state of FDN
          return true;    // Return true if preparation was successful
        }
        // Process a block of samples through the FDN
        bool Process(
          const T* in,                  // Input samples
          T* const outL,                // Left channes output stream
          T* const outR,                // Right channel output stream
          size_t M) noexcept            // Number of samples to process
        {
          for (size_t n=0;n<M;++n)
          {
          T x=T(0);
          if (in)                       // We have an input data stream?
          {                             // Yes, so we are going to
            predelay.Write(in[n]);      // Write the sample to delay line.
            x=predelay.Read();          // Read the sample from the pre-delay line
          }
          const T wet=wetMix.load(std::memory_order_relaxed);
          const T dry=static_cast<T>(1)-wet; // Calculate wet and dry mix
          // Gather current feedback outputs
          for (size_t i=0;i<Ntaps;++i)  // For each filter tap...         
          {                             // Circulate sample through delay-filter graph.
            T r=dls[i].Read();             // Index into sample we are processing.
            T d=dampLP[i].ProcessSample(r); // Damp the sample with the low-pass filter
            lastOut[i]=shelf[i].ProcessSample(d);// Shelf-it
          }
          // Mix through feedback matrix
          
          std:: array<T,Ntaps> feed{};
          for (size_t i=0;i<Ntaps;++i)            // For each delay line
            for (size_t j=0;j<Ntaps;++j)          // For each feedback matrix element
              feed[i]+=fbmtx[i][j]*lastOut[j]; // Mix the outputs through the feedback matrix
          // Add input + feedback into each delay line
          for (size_t i=0;i<Ntaps;++i)            // For each delay line
          {
            dls[i].Write(x+feed[i]); // Write the input + feedback into the delay line
            dls[i].Tick();           // Advance the delay line by one sample
          }
          // Simple stereo tap: even indices -> L, odd indices -> R
          T yL{}, yR{}; // Initialize output samples
          for (size_t i=0;i<Ntaps;++i)
          ((i&1)?yR:yL)+=lastOut[i]; // Add the outputs to left or right channel based on index
          const T norm=static_cast<T>(2)/static_cast<T>(Ntaps); // Normalization factor
          yL*=norm;yR*=norm; // Normalize the output samples
          // Wet/Dry mix
          if (in)
          {
          outL[n]=dry*in[n]+wet*yL;
          outR[n]=dry*in[n]+wet*yR; // Mix the input with the wet output
          }
          else
          {
          outL[n]=wet*yL; // If no input, just output the wet signal
          outR[n]=wet*yR;
          }
        }
          return true; // Return true to indicate successful processing
        }              
        void Clear(void) noexcept
        {
          for (size_t i=0;i<Ntaps;++i)
            shelf[i].Reset();                 // Clear the shelving filter state
          predelay.Clear();           // Clear the pre-delay line
          for (size_t i=0;i<Ntaps;++i)
            dls[i].Clear();                 // Clear the delay line buffer
          for (size_t i=0;i<Ntaps;++i)
            dampLP[i].Reset();                 // Clear the damping filter state
          std::fill(lastOut.begin(), lastOut.end(), T(0)); // Zero the last output samples
        }
        // Set the wet mix level (0.0 to 1.0)
        void SetWetMix(T mix) noexcept
        {
          if (mix<static_cast<T>(0)) mix=static_cast<T>(0);
          else if (mix>static_cast<T>(1)) mix=static_cast<T>(1);
          wetMix.store(mix,std::memory_order_relaxed); // Store the wet mix level
        }
        // Get the current wet mix level
        T GetWetMix(void) const noexcept
        {
          return wetMix.load(std::memory_order_relaxed); // Load the wet mix level
        }
        void SetFeedbackMatrix(detail::MatrixKind k, unsigned seed=0) noexcept
        {
          auto M=detail::MakeMatrix<T,Ntaps>(k,seed);
          constexpr T g=T(0.9);      // Scalar
          for (auto& row:M) 
            for (auto& v: row)           // Find element with max vale and normalize
               v*=g;
          fbmtx=M;                      // Assign normalized matrix.
        }
        
        // Set the feedback matrix (must be orthogonal or unitary)
        void SetFeedbackMatrix(const std::array<std::array<T,Ntaps>,Ntaps>& m) noexcept
        {
          auto M=detail::MakeMatrix<T,Ntaps>(k,seed);
          constexpr T g=T(0.9);      // Scalar
          for (auto& row:M) 
            for (auto& v: row)           // Find element with max vale and normalize
               v*=g;
          fbmtx=M;                      // Ass
          fbmtx=m; // Set the feedback matrix to the provided matrix
        }
        void SetFeedBackMatrix(const sig::spectral::Matrices<T>& M)
        {
          auto M=detail::MakeMatrix<T,Ntaps>(k,seed);
          constexpr T g=T(0.9);      // Scalar
          for (auto& row:M) 
            for (auto& v: row)           // Find element with max vale and normalize
               v*=g;
          fbmtx=M;                      // Ass
          fbmtx=detail::toStdArray<T,Ntaps>(M); // Convert the spectral matrix to std::array
        }
        // Get the current feedback matrix
        const std::array<std::array<T,Ntaps>,Ntaps>& GetFeedbackMatrix(void) const noexcept
        {
          auto M=detail::MakeMatrix<T,Ntaps>(k,seed);
          constexpr T g=T(0.9);      // Scalar
          for (auto& row:M) 
            for (auto& v: row)           // Find element with max vale and normalize
               v*=g;
          fbmtx=M;                      // Ass
          return fbmtx; // Return the current feedback matrix
        }
        // Get the delay lines        
        void SetPreDelaySeconds(double secs) noexcept
        {
            predelay.SetDelay(secs); // Set the pre-delay time in seconds
        }
        
        const std::array<FarrowDelayLine<T>,Ntaps>& GetDelayLines(void) const noexcept
        {
          return dls; // Return the array of delay lines
        }
        void SetDelay(const std::array<double,Ntaps>& secs) noexcept
        {
            for (size_t i=0;i<Ntaps;++i)          // For each delay line
                dls[i].SetDelay(secs[i]); // Set the delay time in seconds
        }
        void SetFractionalDelay(const std::array<double,Ntaps>& secs) noexcept
        {
            for (size_t i=0;i<Ntaps;++i)          // For each delay line
                dls[i].SetFractionalDelay(secs[i]); // Set the fractional delay time in seconds
        }
       void SetDamperCutoffs(const std::array<double,Ntaps>& freq) noexcept
        {
          for (size_t i=0;i<Ntaps;++i)          // For each damping filter
          {
            dampLP[i].Reset();
            dampLP[i].Prepare(fs,freq[i]); // Set the cutoff frequency of the damping filter  
          }
            
        }
        void SetDamperCutoff(size_t idx, double freq) noexcept
        {
          if (idx < 0 || idx >= Ntaps) return; // Sanitize index
          dampLP[idx].Prepare(fs,freq); // Set the cutoff frequency of the specified damping filter
        }
        void SetShelfCutoffs(const std::array<double,Ntaps>& freq) noexcept
        {
            for (size_t i=0;i<Ntaps;++i)          // For each shelving filter
                shelf[i].SetCutoffFrequency(freq[i]); // Set the cutoff frequency of the shelving
        }
        void SetShelfCutoff(int idx, double freq) noexcept
        {
            if (idx < 0 || idx >= Ntaps) return; // Sanitize index
            shelf[idx].SetCutoffFrequency(freq); // Set the cutoff frequency of the specified shelving filter
        }
        void SetShelfSlopes(const std::array<double,Ntaps>& slope) noexcept
        {
            for (size_t i=0;i<Ntaps;++i)          // For each shelving filter
                shelf[i]=filterFactory.LowShelf(fs,shfc,shboost,slope); // Set the slope of the shelving filter
        }
        void SetShelfSlope(int idx, double slope) noexcept
        {
            if (idx < 0 || idx >= Ntaps) return; // Sanitize index
            shelf[idx]=filterFactory.LowShelf(fs,shfc,shboost,slope); // Set the slope of the specified shelving filter
        }
      private:
        double fs{48000.0};
        size_t bs{256};
        FarrowDelayLine<T> predelay;
        double predel{1.3}; // Pre-delay time in seconds
        std::array<FarrowDelayLine<T>,Ntaps> dls;
        std::array<BiQuad<T>, Ntaps> shelf;
        FilterFactory<float> filterFactory;
        std::array<OnePole<T>, Ntaps> dampLP; // Damping filters for each delay line
        double dfc{5000.0};
        double shfc{500.0}; // Shelf cutoff frequency
        double shboost{0.0}; // Shelf boost in dB
        double slope{1.0}; // Shelf slope (1.0 = 6 dB/octave)
        std::array<std::array<T,Ntaps>,Ntaps> fbmtx; // Feedback matrix
        std::array<T,Ntaps> lastOut{}; // Last output samples from each delay line
        std::atomic<T> wetMix{static_cast<T>(0.5)}; //

};
}
