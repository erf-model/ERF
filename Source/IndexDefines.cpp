#include "IndexDefines.H"
#include <AMReX_SPACE.H>

namespace indxmap {

indxmap_values imap;

void
init()
{
  AMREX_D_PICK(imap.upassMap[0] = UMY; imap.qpassMap[0] = QV; imap.upassMap[1] = UMZ;
               imap.qpassMap[1] = QW; int curMapIndx = 2;, imap.upassMap[0] = UMZ;
               imap.qpassMap[0] = QW; int curMapIndx = 1;, int curMapIndx = 0;);
  for (int i = 0; i != NUM_ADV; ++i) {
    imap.upassMap[curMapIndx] = i + UFA;
    imap.qpassMap[curMapIndx] = i + QFA;
    curMapIndx++;
  }
}
} // namespace indxmap
