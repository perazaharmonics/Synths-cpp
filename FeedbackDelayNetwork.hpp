 /*
 * * FeedBackDelayNetwork.hpp
 * * -------------------------------------------------------------
 * * A high-quality, compile-time-sized FDN core using FarrowDelayLine
 * * for smooth fractional delays and modulation.  Designed as the
 * * foundation for echo, reverb, chorus/flanger, tremolo, etc.
 * * Before reading this modules, please read the DelayBranch.hpp
 * * as that single waveguide branch  is the basis for the FDNs.
 * * 
 * * NOTE: FDNs are not only useful for reverbs and echoes, but
 * *       because they are a structure of interconnected waveguides,
 * *       they can be used for a variety of applications concerning
 * *       physical modeling, spatialization, and even machine learning.
 * *       Thus they are seminal to the field of DSP, simulation software
 * *       like CADs and HFSS, RF and Wireless Communications.
 * *
 * * Template parameters
 * * -------------------
 * *   T         : sample type (float / double)
 * *   MaxLen    : maximum length of the delay line in samples
 * *   K         : Farrow filter order (default 3)
 * *   P         : Thiran filter order (default 3)
 * *   Ntaps     : number of delay branches (tap-off points) in the FDN (default 8)
 * * Key Features
 * * ------------
 * * -> Fixed Branch for every line
 * * -> Arbitrary orthogonal / unitary feedback matrix
 * * -> Per-line one-pole damping filters (high-shelf and/or LPF)
 * * -> Smooth parameter changes (thread-safe)
 * * -> Wet / Dry mix, Output tap-matrix for stereo or multichannel
 * *
 * * Uses:
 * *  -> Resonator Bank && Physical Modeling of Waveguides:
 * *       Choose the tap-lengths and feedback matrix gains to mimic modal frequencies,
 * *       an the FDN becomes a resonator bank.
 * *
 * * ->  Decorrelation and Spatialization: 
 * *        Use Hadamard or Householder matrices to decorrelate the taps. Feed the same dry signal
 * *        at different delay lengths (phases) to the FDN (waveguide network). The output
 * *        will be decorrelated of "early-reflection" cues. Good for virtual surround sound, headphone spatialization
 * *        and widening mono sources like vocals.
 * *
 * * -> Multi-Tap Comb Filter Network: 
 * *      Collapse or zero-out some of the waveguide branches to create a "comb-filter forest" structure.
 * *      Good for phasing signals (think beamforming) and adaptive resonant filtering.
 * * 
 * * -> Reservoir Computing/Echo-State Networks:
 * *     Use the FDN as a fixed-recurrent network "reservoir" and learn a simple linear readout on top.
 * *     This allows for real-time pattern recognition and classification, simple predictive modelling, or
 * *     rudimentary Machine Learning tasks using the FDN as a dynamic feature extractor.
 * *
 * * -> Multipath Channel Simulation:
 * *     Because the FDN is just a structure of interconnected Waveguides, the FDN is used for
 * *     communication-channel multipath echo and fades simulation (Underwater Acoustics, RF, MMwave and Microwave channels)
 * *     for testing equalizers or channel-estimation algorithms.
 * *
 * * -> Noise Shaping and Dithering:
 * *     Use the FDN with a carefully designed feedback matrix to make the FDN act as a high-order multi-tap
 * *     noise shaper that can push quantization noise and quantization error to the signal's null-space (i,e, the noise floor).
 * *
 * * -> Multi-Band All-Pass Waveguide Network:
 * *      Use the FDN in a cascade with another FDN tuned to different frequency bands. This will allow you to create
 * *      Maximally Flat All-Pass Waveguide resonant to different frequency bands, and perform controlled phase-altering crossovers
 * *      by tuning the FDN, allowing for mid/side manipulation, dynamic phase-alignment, transient shaping inserts,
 * *      and beamforming.
 * * 
 * * Author:
 * * JEP J. Enrique Peraza
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
#include "Matrices.hpp"      // Matrix base class 
#include "MatrixSolver.hpp"  // This is just a really powerful eigenvector/eigenvalue solver. (If you know why we need it, you know....)

namespace sig::wg {
   
  namespace detail                        // Helper namespace 
  {
    template <typename T, size_t N> // N filter taps.
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

    template <typename T, size_t N>
    constexpr auto hadamard() noexcept
    {
        static_assert((N & (N - 1)) == 0, "Hadamard requires power-of-two N");
        std::array<std::array<T, N>, N> m{};
        m[0][0] = static_cast<T>(1);
        for (size_t size = 1; size < N; size <<= 1)
          for (size_t i = 0; i < size; ++i)
            for (size_t j = 0; j < size; ++j)
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
    template <typename T, size_t N>
    inline auto randomOrthogonal(unsigned seed = 0xCEFAFD) noexcept
    {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<T> dist;
        std::array<std::array<T, N>, N> a{};
        for (auto& row:a)
            for (auto& v:row) v=dist(gen);
        // Gram–Schmidt QR ? Q
        for (size_t k=0;k<N;++k)
        {
          for (size_t i=0;i<k;++i)
          {
            T dot{};
            for (size_t n=0;n<N;++n) dot+=a[k][n]*a[i][n];
            for (size_t n=0;n<N;++n) a[k][n]-=dot*a[i][n];
          }
          T norm{};
          for (auto v:a[k]) norm+=v*v;
          norm = std::sqrt(norm);
          for (auto& v:a[k]) v/=norm;
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

    template<typename T, size_t N>
    constexpr std::array<std::array<T,N>,N> toStdArray(const Matrices<T>& M)
    {
        if (M.Rows() != N || M.Cols() != N)
            throw std::invalid_argument{"Matrix size ? FDN order"};

        std::array<std::array<T,N>,N> out{};
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j)
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
          fbmtx=sig::wg::detail::householder<T, Ntaps>(); // Init feedback matrix
          mtxk=sig::wg::detail::MatrixKind::Householder; // Track kind
          for (size_t i=0;i<Ntaps;++i)
              dls[i].Prepare(static_cast<size_t>(defaultDelays[i]), /*Thiran*/0.0f,/*Farrow*/0.0f); // Use Prepare instead of SetDelay
      }
      ~FeedbackDelayNetwork(void) noexcept = default;

      // Prepare the FDN with a given delay time and damping factor
      bool Prepare(double fs,size_t bs)
      {                                 // Prepare the FDN with a given sample rate and block size
        if (fs<=0.0||bs<1) return false; // If sample
        this->fs=fs;                    // Set the sample rate
        this->bs=bs;                    // Set the block size
        if (predel>0.0)                 // Cealr predelay state
        {
          predelay.Prepare(static_cast<size_t>(predel),0.0f,0.0f); // Prepare the predelay
          predelay.Clear();
        }// Clear FDN state.
        for (auto& dl:dls)             // For each delay line
          dl.Clear();                  // Clear the delay line state
        for (size_t i=0;i<Ntaps;++i)   // For each tap (matrix row)
        {                              // Check if we have fractional delays set.
          if (mut[i]>0.0f)             // Any Thiran fraction requested?
            dls[i].SetFractionalDelay(mut[i]); // Set Thiran fractional delay
          if (muf[i]>0.0f)             // Any Farrow modulation requested?
            dls[i].SetMuFarrow(muf[i]); // Set fractional delay for Farrow
        }
        for (size_t i=0;i<Ntaps;++i)    // For each tap (matrix row)
          shelf[i]=filterFactory.LowShelf(fs,shfc,shboost,slope);
        for (size_t i=0;i<Ntaps;++i)    // For each tap, 
        {                               // Set damper filter and prepare DelayBranch
          size_t idelay=dls[i].GetDelay(); // Get integer delay in samples from DelayBranch
          dampLP[i] = filterFactory.OnePoleLP(fs, dfc); // This filter from the factory.
          dampLP[i].Prepare(fs, bs);    // Prepare the damper filter
          dls[i].Prepare(idelay,mut[i],muf[i]); // Prepare the DelayBranch
        }                                // Done processing through 
        size_t maxlat=0;                // Maximum latency in samples
        for (auto& d:dls)                // For each delay line
        {
          size_t D=d.GetDelay();       // Get the group delay in samples
          int G=d.GroupDelay(D, P, K); // Get the group delay in samples
          maxlat=std::max(maxlat,D+G); // Update maximum latency
          // ------------------------- //
          // Prime the whole FDN by running 0 samples for length of group delay
          // ------------------------- //
          T dummyL,dummyR;             // Dummy output for the delay line
          T* dummyIn=nullptr;          // Dummy input buffer
          dummyIn=new T[maxlat];       // Create dummy input buffer
          for (size_t i=0;i<maxlat;++i)// For the group delay length...
          {                            // Pump zeroes to prime the FDN
            dummyIn[i]=T(0); // Fill the dummy input buffer with zeros
            Process(&dummyIn[i], &dummyL, &dummyR, 1); // Process zeroes down the FDN
          }                             // Done priming the FDN.
          delete[] dummyIn;             // Delete the dummy input buffer
          dummyIn=nullptr;              // Clear the dummy input buffer
        }                               // Done priming the FDN.
        return true;                    // Return true if preparation was successful
      }                                 // Prepare the FDN with a given delay time and damping factor
      // Process a block of samples through the FDN
      bool Process(
        const T* in,                    // Input buffer (nullptr for dry output)
        T* const outL,                  // Left output buffer pointer to stream
        T* const outR,                  // Right output buffer pointer to stream
        size_t M) noexcept              // Process a block of samples through the FDN 
      {                                 // ----------- Process --------------------
        for (size_t n=0; n<M;++n)       // For each sample in the block
        {
        bool nofb=true;             // Record if feedback state is empty
        for (size_t i=0;i<Ntaps&&nofb;++i)// For each tap (line)
          for (size_t j=0;j<Ntaps;++j) // For each tap (column)
            if (fbmtx[i][j]!=T(0)) nofb=false; // For the feedback matrix....
        T x=T(0);                   // Initialize output sample
        if (in!=nullptr)            // Is the input null?
        {                           // No, we can go on
          if (wetMix.load(std::memory_order_relaxed)>=T(1.0)&&nofb)
            x=in[n];                // Copy the input sample
          if (predel>0.0)           // User wants predelay?
          {
            predelay.Write(in[n]); // Write the input sample to predelay
            predelay.Propagate(1); // Propagate the predelay
            x=predelay.Read();     // Read the predelay output
          }
          else                     // No predelay, just copy input
            x=in[n];               // Copy the input sample
        }
        const T wet=wetMix.load(std::memory_order_relaxed);
        const T dry=static_cast<T>(1)-wet;
        bool nodelay=true;       // True if no delay 
        for (auto& d:dls)
          if (d.GetDelay()>0) {nodelay=false;break;}
          // If no feedback and wet=1.0, just copy input to output
        if (in&&wet==T(1.0)&&nofb&&nodelay&&predel<=T(0.0))
        {                             // ~~~~~~~~~~~~~~~~~~~~~~
          outL[n]=x;                  // Laaaaaaazy
          outR[n]=x;                  //          Coooooooopy
        }                            // Done doing lazy copy.
        for (size_t i=0;i<Ntaps;++i)  // For the number of coeffs in filter
        {                             // Feed into the filter blocks
          // 1. Read the delay line output
          T r=dls[i].Read();         // Read the delay line output
          // ----------------------- //
          // Apply shelf and damping filters
          // ----------------------- //
          T d=(dampfc[i]>0.0&&dampfc[i]<fs*0.5)?dampLP[i].ProcessSample(r):r;
          // Apply the shelf and damper only if below Nyquist limit!
          lastOut[i]=(shelffc[i]>0.0&&shelffc[i]<fs*0.5)?shelf[i].ProcessSample(d):d;
        }                             // Done applying filtersssz.
        // -------------------------- //
        // 2. Mix through feedback matrix
        // -------------------------- //
        std::array<T,Ntaps> feed{};  // Initialize feedback array
        for (size_t i=0;i<Ntaps;++i)
        {
          T acc=0;
          for (size_t j=0;j<Ntaps;++j)
            acc+=fbmtx[i][j]*lastOut[j]; // Mixed output
          feed[i]=acc;                  // Store the mixed output
        }
        // -------------------------- //
        // 3. Data dumping into each delay line
        // -------------------------- //
        /// 3. Write the input + feedback to each delay line + propagate again.
        for (size_t i=0;i<Ntaps;++i)
        {
          dls[i].Write(x+feed[i]);    // Write the input + feedback to the delay line
          dls[i].Propagate(1);         // Propagate the delay line
        }
        // -------------------------- //
        // 4. Simple stereo tap: even ? L, odd ? R
        // -------------------------- //
        T yL=0,yR=0;                  // Initialize left and right outputs
        for (size_t i=0;i<Ntaps;++i)  // Corculate through blth channels
        {
          // tap position 0..1
          T t=Ntaps>1?T(i)/T(Ntaps-1):T(0);
          // Map to angle
          T th=t*static_cast<T>(M_PI_2); // Map to angle
          // constant power gain
          T gL=std::cos(th); // Left gain
          T gR=std::sin(th); // Right gain
          // Apply gain to the output
          yL+=lastOut[i]*gL; // Left output
          yR+=lastOut[i]*gR; // Right output
        }
        // Normalize energy across N taps
        T invG=T(1)/std::sqrt(static_cast<T>(Ntaps));
        yL*=invG; // Normalize left output
        yR*=invG; // Normalize right output
        if (in!=nullptr)
        {
          outL[n]=dry*in[n]+wet*yL;
          outR[n]=dry*in[n]+wet*yR; // Mix dry and wet outputs
        }
        else
        {
          outL[n]=wet*yL;           // Write wet output to left channel
          outR[n]=wet*yR;           // Write wet output to right channel
        }
        }
       return true;                   // Return true if processing was successful
      }                               // Process a block of samples through the FDN
      void Clear(void) noexcept       // Reset the FDN state
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
          fbmtx=detail::MakeMatrix<T,Ntaps>(k,seed);
          mtxk=k;
      }
      void SetFeedbackMatrix(const std::array<std::array<T,Ntaps>,Ntaps>& m) noexcept
      {
          fbmtx=m;
          mtxk=detail::MatrixKind::Identity;
      }
      const std::array<std::array<T,Ntaps>,Ntaps>& GetFeedbackMatrix() const noexcept { return fbmtx; }

      using Matrices = std::vector<std::vector<T>>;
      void SetFeedBackMatrix(const Matrices m) { fbmtx = detail::toStdArray<T,Ntaps>(m); }

      // Delay helpers
      void SetPreDelay(double samps) noexcept
      {
          predel = samps;
          predelay.Prepare(static_cast<size_t>(predel),0.0f,0.0f); // Use Prepare instead of SetDelay
      }
      // --------------------------------------------------------------------- //
      //  Compatibility shims – keep old client code alive with the new engine //
      // --------------------------------------------------------------------- //
      inline void SetDelay(const std::array<double,Ntaps>& delays) noexcept
      {
          defaultDelays=delays; // Set the default delays
          for (size_t i=0;i<Ntaps;++i)           // For each delay line
              dls.Prepare(static_cast<size_t>(std::floor(delays[i])),/*mu Thiran*/0.0f,/*mu Farrow*/0.0f); // Use Prepare to set delay and fractional delay
      }
      inline void Tick(void) noexcept
      {
        for (auto& dl:dls)      // For # of delay branches in system
          dl.Propagate(1);        // Make it Tick().....
      }
      const std::array<Branch,Ntaps>& GetDelayLines(void) const noexcept { return dls; }
      void SetFractionalDelay(const std::array<double,Ntaps>& s) noexcept
      {                               // Set Thiran fractional delays
        for (size_t i=0;i<Ntaps;++i)  // For the number of taps....
        {                             // Tune Thiran filters
          mut[i]=static_cast<double>(s[i]);// Save the user's values.
          dls[i].SetFractionalDelay(s[i]);// Set the Thiran fractional delay
        }                             // Done tuning
      }                               // Set Thiran fractional delays
      void SetMuFarrow(const std::array<double,Ntaps>& s) noexcept
      {                                // Set Farrow fractional delays.
        for (size_t i=0;i<Ntaps;++i)   // For the number of taps...
        {                              // Tune Farrow modulation filters
          muf[i]=static_cast<double>(s[i]);// Save the user's values.
          dls[i].SetMuFarrow(s[i]);    // Set the Farrow fractional delay
        }                              // Done tuning
      }                                // Set Farrow fractional delays.
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
      { 
        defaultDelays=delays; // Set the default delays
        for (size_t i=0;i<Ntaps;++i)
          dls[i].Prepare(static_cast<size_t>(delays[i]), 0.0f, 0.0f); 
      }

  private:
      double fs{48000.0};
      size_t bs{256};
      Branch predelay;
      double predel{0.0};        // Pre-delay time in seconds
      std::array<Branch,Ntaps> dls;
      std::array<double,Ntaps> defaultDelays{}; // Initial delay times in seconds
      std::array<double,Ntaps> mut{};           // Thiran fractional delays
      std::array<double,Ntaps> muf{};           // Farrow fractional delays
      std::array<BiQuad<T>,Ntaps> shelf;
      FilterFactory<float> filterFactory;
      std::array<OnePole<T>,Ntaps> dampLP;
      std::array<double,Ntaps>  dampfc{};
      std::array<double,Ntaps> shelffc{};
      double dfc{fs*0.5};
      double shfc{5000.0};
      double shboost{0.0f};
      double slope{0.0};
      std::array<std::array<T,Ntaps>,Ntaps> fbmtx;
      std::array<T,Ntaps>     lastOut{};
      std::atomic<T>          wetMix{static_cast<T>(0.0)};
      detail::MatrixKind      mtxk = detail::Identity;

      static void Normalize(std::array<std::array<T,Ntaps>,Ntaps>& M) noexcept
      {
          T maxv=T(0);
          for (const auto& row:M) for (const auto& v:row) maxv = std::max(maxv,std::fabs(v));
          if (maxv!=T(0)) for (auto& row:M) for (auto& v:row) v/=maxv;
      }
  };
} // namespace sig::wg
