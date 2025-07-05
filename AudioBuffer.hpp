/**
 * Filename: AudioBuffer.hpp
 * Description: Simple planar (non-interleaved) audio buffer.
 * 
 * Author:
 * JEP J. Enrique Peraza
 */
#pragma once
#include <vector>
#include <cstddef>
#include <cstring>
namespace sig{
template<typename T>
class AudioBuffer{
public:
  // Contructor:
  AudioBuffer(void)=default;
  // Resize; returns false on OOM and leaves old data intact.
  bool Resize(std::size_t nFrames,std::size_t nChannels) noexcept
  {                                     // ---------- Resize --------------- //
      this->frames=nFrames;             // Set the number of frames.
      this->channels=nChannels;         // Set the number of channels.
      this->data.resize(nChannels*nFrames);// Set data vector size.
    return true;                        // Return true if we got here.
  }                                     // ---------- Resize --------------- //
  // Getters:
  [[nodiscard]] std::size_t Frames(void) const noexcept { return frames; }
  [[nodiscard]] std::size_t Channels(void) const noexcept { return channels; }
  // Pointer to channel i (0-based, planar)
  [[nodiscard]] T* Channel(std::size_t i) noexcept {return data.data()+i*frames;}
  [[nodiscard]] const T* Channel(std::size_t i) const noexcept {return data.data()+i*frames;}

  // Clear the buffer
  void Clear(void) noexcept {std::memset(data.data(),0,data.size()*sizeof(T));}
 private:
   std::vector<T> data;                 // The audio signal.
   std::size_t frames;                  // The number of frames.
   std::size_t channels;                // The number of channels. 
};                                      // Class: AudioBuffer
}                                       // namespace sig                    