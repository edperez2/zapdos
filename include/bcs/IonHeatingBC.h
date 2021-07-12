//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADIntegratedBC.h"

class IonHeatingBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  IonHeatingBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _r_units;
  const Real _r;
  const Real _r_ion;
  const Real _ionization_energy;
  const MaterialProperty<Real> & _kb;

  // Coupled variables

  const ADVariableGradient & _grad_potential;
  const ADVariableValue & _em;
  const ADVariableValue & _potential;
  const ADVariableValue & _mean_en;
  std::vector<const ADVariableValue *> _ip;
  std::vector<const ADVariableGradient *> _grad_ip;

  const ADMaterialProperty<Real> & _muem;
  const MaterialProperty<Real> & _massem;
  const MaterialProperty<Real> & _e;
  std::vector<const MaterialProperty<Real> *> _sgnip;
  std::vector<const ADMaterialProperty<Real> *> _muip;
  std::vector<const ADMaterialProperty<Real> *> _Tip;
  std::vector<const MaterialProperty<Real> *> _massip;
  const MaterialProperty<Real> & _se_coeff;
  const MaterialProperty<Real> & _se_energy;
  const ADMaterialProperty<Real> & _mumean_en;

  Real _a;
  Real _b;
  ADReal _v_thermal;
  ADReal _ion_flux;
  ADReal _n_gamma;
  ADReal _energy_flux;

  unsigned int _num_ions;
};
