#pragma once
#include <cstddef>   // std::size_t

namespace sig::wg
{
// ------------------------------------------------------------
// A Sample is just a strongly typed scalar  alias to T
// ------------------------------------------------------------
template<typename T> using Sample = T;

// ------------------------------------------------------------
// Base class for every processing element
// ------------------------------------------------------------
struct Node
{
    virtual ~Node() = default;
    virtual void Prepare(std::size_t /*nFrames*/) noexcept {}   // optional
    virtual void Propagate(std::size_t nFrames)   noexcept = 0; // must implement
};

// ------------------------------------------------------------
// Forward declaration for DelayBranch (real template defaults in its definition)
template<
    typename T,
    std::size_t MaxLen,
    std::size_t FarOrder,
    std::size_t ThOrder
 >
 class DelayBranch;
} // namespace sig::wg
