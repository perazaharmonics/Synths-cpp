/* 
 * * 
 * *
 * * Filename: InstrumentGraph.hpp
 * * Description:
 * * This file contains a InstrumentGraph class that stores the connection
 * * pf the Node* object in topological order.
 * * It is used to model the acoustic properties of frequency travelling through
 * * a networks of waveguide and scattering elements. I.E., Nodes with a particular
 * * property that performs a computation on the input sample and produces an output sample, which is then
 * *  deinterpolated and projected unto the return delay line of the waveguide.
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <vector>
#include "WGTypes.hpp"
namespace sig::wg
{
  class InstrumentGraph
  {
  public:
    void AddNode(Node* n) noexcept { this->nodes.push_back(n); } // Add a node to the graph.
    void Propagate(size_t nFrames) noexcept
    {                                   // ----------- Propagate ----------------- //
      for (auto* n:this->nodes)         // For every node in the graph...
        n->Propagate(nFrames);          // Propagate the samples through the node.
    }                                   // ----------- Propagate ----------------- //
    inline const std::vector<Node*>& GetNodes(void) const noexcept { return this->nodes; } // Get the nodes in the graph.
    inline void Clear(void) noexcept { this->nodes.clear(); } // Clear the nodes in the graph.
    inline size_t Size(void) const noexcept { return this->nodes.size(); } // Get the number of nodes in the graph.
  private:
    std::vector<Node*> nodes; // Vector to store the nodes in the graph.
  // This vector stores pointers to Node objects, which are the elements of the graph.
  // Each Node object represents a component in the waveguide network, such as a string element or a cylindrical bore.
  // The nodes are stored in topological order, meaning that they can be processed in a way that respects their dependencies.
  // The Propagate method iterates through each node and calls its Propagate method, allowing the waveguide network to process samples in a sequential manner.
  // This design allows for flexible and dynamic construction of waveguide networks, where nodes can be added or removed as needed.
  };
}