// Created by Pankaj Jha on 6/1/21.

#include "ERF.H"
using namespace amrex;

enum class BCType_User{
  FreeSlip, NoSlip, AnotherBC
};

// Dirichlet BC will need as an argument the fixed value of the variable at the boundary
Real DirichletBC (const int &i, const int &j, const int &k,
                        const Array4<Real>& var, const Real &boundaryVal);

// Neumann BC will need as an argument the gradient of the variable at the boundary
Real NeumannBC (const int &i, const int &j, const int &k,
                  const Array4<Real>& var, const Real &gradientVal);

// Neumann BC will need as an argument the relevant value of the variable at the boundary
Real AnotherBC (const int &i, const int &j, const int &k,
                const Array4<Real>& var, const Real &relevantVal);



void applyBC (const int &i, const int &j, const int &k, const enum BCType_User &bcUser)
{
  // Declared it here just for demo. The BC should be applied on the actual Array4<Real> corresponding to a multi-fab
  Array4<Real> u, v, w, theta, T, p; // Possibly pass u, v etc. to this function

  switch (bcUser) {
    case BCType_User::FreeSlip:
      Real u_bdry, v_bdry, T_bdry, p_relevantVal; // These should be available from input
      DirichletBC(i, j, k, u, u_bdry);
      DirichletBC(i, j, k, v, v_bdry);
      DirichletBC(i, j, k, T, T_bdry);
      AnotherBC(i, j, k, p, p_relevantVal);
      break;
    case BCType_User::NoSlip:
      Real u_grad, v_grad, T_grad; // These should be available from input
      NeumannBC(i, j, k, u, u_grad);
      NeumannBC(i, j, k, v, v_grad);
      NeumannBC(i, j, k, T, T_grad);
    case BCType_User::AnotherBC:
      Real w_someVal, theta_someVal; // These should be available from input
      AnotherBC(i, j, k, w, w_someVal);
      AnotherBC(i, j, k, theta, theta_someVal);
      break;
    default:
      amrex::Abort("Error: Boundary condition is unrecognized");
  }
}