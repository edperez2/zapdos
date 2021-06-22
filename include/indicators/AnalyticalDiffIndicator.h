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

#include "ElementIntegralIndicator.h"

class AnalyticalDiffIndicator : public ElementIntegralIndicator
{
public:
  AnalyticalDiffIndicator(const InputParameters & parameters);

  static InputParameters validParams();

  virtual ~AnalyticalDiffIndicator(){};

protected:
  virtual Real computeQpIntegral();

  const Function & _func;
};
