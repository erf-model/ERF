#include <TimeIntegration.H>

using namespace amrex;

int
ComputeGhostCells(const int& spatial_order) {
  int nGhostCells;

  switch (spatial_order) {
    case 2:
      nGhostCells = 1;
      break;
    case 3:
      nGhostCells = 2;
      break;
    case 4:
      nGhostCells = 2;
      break;
    case 5:
      nGhostCells = 2;
      break;
    case 6:
      nGhostCells = 3;
      break;
    default:
      amrex::Error("Must specify spatial order to be 2,4, or 6");
  }

  return nGhostCells;
}
