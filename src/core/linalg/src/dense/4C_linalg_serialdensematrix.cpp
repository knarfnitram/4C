// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  B = alpha*A + beta*B                                                |
 *----------------------------------------------------------------------*/
void Core::LinAlg::update(double alpha, const Core::LinAlg::SerialDenseMatrix& A, double beta,
    Core::LinAlg::SerialDenseMatrix& B)
{
  B.scale(beta);
  Core::LinAlg::SerialDenseMatrix Acopy(A);
  Acopy.scale(alpha);
  B += Acopy;
}

/*----------------------------------------------------------------------*
 |                                                         vlf 06/07    |
 | recursive computation of determinant of a  matrix using Sarrus rule  |
 | (uses long double to boost accuracy). Do not use for n > 4.          |
 *----------------------------------------------------------------------*/
long double Core::LinAlg::det_long(const Core::LinAlg::SerialDenseMatrix& matrix)
{
  if (matrix.numCols() == 1)
  {
    return matrix(0, 0);
  }
  else if (matrix.numCols() == 2)
  {
    long double out_det;
    out_det = ((long double)(matrix(0, 0)) * (long double)(matrix(1, 1))) -
              ((long double)(matrix(0, 1)) * (long double)(matrix(1, 0)));
    return out_det;
  }
  else if (matrix.numCols() > 2)
  {
    long double out_det = 0;
    int sign = 1;
    for (int i_col = 0; i_col < matrix.numCols(); i_col++)
    {
      SerialDenseMatrix temp_matrix(matrix.numCols() - 1, matrix.numCols() - 1);
      for (int c_col = 0; c_col < i_col; c_col++)
      {
        for (int row = 1; row < matrix.numCols(); row++)
          temp_matrix(row - 1, c_col) = matrix(row, c_col);
      }
      for (int c_col = i_col + 1; c_col < matrix.numCols(); c_col++)
      {
        for (int row = 1; row < matrix.numCols(); row++)
          temp_matrix(row - 1, c_col - 1) = matrix(row, c_col);
      }
      out_det = out_det + ((long double)(sign) * (long double)(matrix(0, i_col)) *
                              (long double)(det_long(temp_matrix)));
      sign *= -1;
    }
    return out_det;
  }
  else
    return 0;
}

/*----------------------------------------------------------------------*
 |   Set matrix components to zero                                      |
 |   this(0:Length) = 0.0                                               |
 *----------------------------------------------------------------------*/
void Core::LinAlg::zero(Core::LinAlg::SerialDenseMatrix& mat, int Length)
{
  int cnt = 0;
  for (int j = 0; j < mat.numCols(); ++j)
  {
    for (int i = 0; i < mat.numRows(); ++i)
    {
      if (cnt++ < Length)
        mat(i, j) = 0.0;
      else
        break;
    }
  }
}

/*----------------------------------------------------------------------*
 |   Fill the SerialDenseMatrix column-wise with vector data            |
 *----------------------------------------------------------------------*/
void Core::LinAlg::copy(const double* vec, Core::LinAlg::SerialDenseMatrix::Base& mat)
{
  int cnt = 0;
  for (int j = 0; j < mat.numCols(); ++j)
  {
    for (int i = 0; i < mat.numRows(); ++i)
    {
      mat(i, j) = vec[cnt++];
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
