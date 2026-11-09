#pragma once
// Setup report of the CRP PES (SystemID CRPlibrary): the data set and how the grid relates to its
// range. The Fortran prints only the library's own messages.

#include <ostream>

#include "potential/crp_potential.hpp"
#include "report/format.hpp"

namespace spodvr {

// The data set (symmetry, cuts, ranges, switching, energy shift, warnings) and the grid points
// that are clamped or have no data (no Z points: analysis-only runs, no 6D potential).
void print_crp_data_set(std::ostream& out, const Style& style, const CrpDataSetSummary& summary,
                        const CrpRangeReport& range);

}  // namespace spodvr
