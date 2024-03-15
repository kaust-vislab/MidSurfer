// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
#include "vtkRungeKutta2Linefield.h"

#include "vtkFunctionSet.h"
#include "vtkObjectFactory.h"
#include <vtkLogger.h>
#include <vtkMath.h>

vtkStandardNewMacro(vtkRungeKutta2Linefield);

vtkRungeKutta2Linefield::vtkRungeKutta2Linefield() = default;

vtkRungeKutta2Linefield::~vtkRungeKutta2Linefield() = default;

void vtkRungeKutta2Linefield::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

// Calculate next time step
int vtkRungeKutta2Linefield::ComputeNextStep(double *xprev, double *dxprev, double *xnext, double t,
                                             double &delT, double &delTActual, double, double, double, double &error, void *userData)
{
  int i, numDerivs, numVals;

  delTActual = 0.;
  error = 0.0;

  if (!this->FunctionSet)
  {
    vtkErrorMacro("No derivative functions are provided!");
    return NOT_INITIALIZED;
  }

  if (!this->Initialized)
  {
    vtkErrorMacro("Integrator not initialized!");
    return NOT_INITIALIZED;
  }

  numDerivs = this->FunctionSet->GetNumberOfFunctions();
  numVals = numDerivs + 1;
  for (i = 0; i < numVals - 1; i++)
  {
    this->Vals[i] = xprev[i];
  }
  this->Vals[numVals - 1] = t;

  // Obtain the derivatives dx_i at x_i
  if (dxprev)
  {
    for (i = 0; i < numDerivs; i++)
    {
      this->Derivs[i] = dxprev[i];
    }
  }
  else if (!this->FunctionSet->FunctionValues(this->Vals, this->Derivs, userData))
  {
    memcpy(xnext, this->Vals, (numVals - 1) * sizeof(double));
    return OUT_OF_DOMAIN;
  }

  // step backward
  double tmpvals[4];
  for (i = 0; i < numVals - 1; i++)
  {
    tmpvals[i] = xprev[i] - this->direction * delT * this->Derivs[i];
  }
  tmpvals[numVals - 1] = t - delT;

  // Obtain the derivatives at x_i - dt * dx_i
  double dxtmp[4];

  if (this->FunctionSet->FunctionValues(tmpvals, dxtmp, userData))
  {
    vtkLog(INFO, "Derivs: " << this->Derivs[0] << ", " << this->Derivs[1]);
    vtkLog(INFO, "dxtmp: " << dxtmp[0] << ", " << dxtmp[1]);
    
    // if (vtkMath::Dot(this->Derivs, dxtmp) < 0.0)
    if ((this->Derivs[0] * dxtmp[0] + this->Derivs[1] * dxtmp[1]) < 0.0)
    {
      vtkLog(INFO, "Reversing direction.");
      this->direction *= -1.0;
    }
  }

  // Calculate x_i using improved values of derivatives
  for (i = 0; i < numDerivs; i++)
  {
    xnext[i] = xprev[i] + this->direction * delT * this->Derivs[i];
  }

  delTActual = delT;

  return 0;
}
