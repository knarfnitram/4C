// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_GRAPH_HPP
#define FOUR_C_LINALG_GRAPH_HPP


#include "4C_config.hpp"

#include <Epetra_CrsGraph.h>
#include <Epetra_Export.h>

#include <memory>

// Do not lint the file for identifier names, since the naming of the Wrapper functions follow the
// naming of the Epetra_Vector

// NOLINTBEGIN(readability-identifier-naming)

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{

  class Graph
  {
   public:
    /// Copy constructor from Epetra_CrsGraph to vector
    explicit Graph(const Epetra_CrsGraph& Source);

    //! Creates a Epetra_CrsGraph object and allocates storage.
    Graph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, const int* NumIndicesPerRow,
        bool StaticProfile = false);

    Graph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int NumIndicesPerRow,
        bool StaticProfile = false);

    // Rule of five: We currently need to take care to make a deep copy of the Epetra_CrsGraph.
    Graph(const Graph& other);

    Graph& operator=(const Graph& other);

    ~Graph() = default;

    //! (Implicit) conversions: they all return references or RCPs, never copies
    const Epetra_CrsGraph& get_ref_of_Epetra_CrsGraph() const { return *graph_; }

    Epetra_CrsGraph& get_ref_of_Epetra_CrsGraph() { return *graph_; }

    std::shared_ptr<Epetra_CrsGraph> get_ptr_of_Epetra_CrsGraph() { return graph_; }


    // Functions

    //! Returns the Column Map associated with this graph.
    const Epetra_BlockMap& ColMap() const { return (graph_->ColMap()); }

    const Epetra_Export* Exporter() { return graph_->Exporter(); };

    int Export(const Epetra_SrcDistObject& A, const Epetra_Export& Exporter,
        Epetra_CombineMode CombineMode, const Epetra_OffsetIndex* Indexor = 0)
    {
      return graph_->Export(A, Exporter, CombineMode, Indexor);
    }

    int FillComplete(void);

    //! Enter a list of elements in a specified global row of the graph.
    int InsertGlobalIndices(int GlobalRow, int NumIndices, int* Indices);

    int InsertGlobalIndices(long long GlobalRow, int NumIndices, long long* Indices);

    //! Make consecutive row index sections contiguous, minimize internal storage used for
    //! constructing graph
    int OptimizeStorage() { return graph_->OptimizeStorage(); }

    //! Remove a list of elements from a specified global row of the graph.
    int RemoveGlobalIndices(int GlobalRow, int NumIndices, int* Indices);

    int RemoveGlobalIndices(long long GlobalRow, int NumIndices, long long* Indices);

    const Epetra_BlockMap& RowMap() const { return graph_->RowMap(); }

    explicit Graph(const std::shared_ptr<Epetra_CrsGraph>& graph);


   private:
    //! The actual Epetra_CrsGraph object.
    std::shared_ptr<Epetra_CrsGraph> graph_;
  };
}  // namespace Core::LinAlg
FOUR_C_NAMESPACE_CLOSE

// NOLINTEND(readability-identifier-naming)

#endif