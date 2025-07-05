/**
 * Filename: IProcessor.hpp
 * 
 * Description: Generic sig block interface (no exceptions, RT-safe)
 * 
 * Author:
 * JEP J. Enrique Peraza
 */
#pragma once
#include <cstddef>
#include "AudioBuffer.hpp"
namespace sig{
template<typename T>
class IProcessor{
public:
  virtual ~IProcessor(void)=default;
  // Called once before the first call to proces().
  // blockSiz: Host I/O buffer size in frames.
  // numCh: number of channels that will appear in every AudioBuffer.
  // Return false on any fatal error (e.g. cannot allocate internal buffers)
  virtual bool Prepare(
    double fs,                          // The sampling frequency.
    std::size_t blockSize,              // The host I/O buf size in frames.
    std::size_t numCh) noexcept=0;      // Number of channels per AudioBuffer.

  // RT-safe processing - **must not allocate** or throw.
  virtual void Process(AudioBuffer<T>& io,std::size_t nFrames) noexcept=0;
  // Reset internal state (called on transport jumps).
  virtual void Reset(void) noexcept {}
};
} // namespace sig.
