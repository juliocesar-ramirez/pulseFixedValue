/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "direction.H"
#include "pulseFixedValueFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include <cmath>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::pulseFixedValueFvPatchField<Type>::t() const
{
    return this->db().time().timeOutputValue();
}

template<class Type>
Foam::scalar Foam::pulseFixedValueFvPatchField<Type>::pulseFraction() const
{
    scalar duracionFraction=duracion_/periodo_;
    scalar cycleFraction=fmod(t()/periodo_,1);
    if (cycleFraction>duracionFraction) {
        return 0.0;;
    }
    else {
        return 1.0;
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Type>
Foam::pulseFixedValueFvPatchField<Type>::pulseFixedValueFvPatchField(
    const fvPatch &p, const DimensionedField<Type, volMesh> &iF)
    : fixedValueFvPatchField<Type>(p, iF), periodo_(0.0),duracion_(0.0), baseValue_(Zero),
      pulseValue_(p.size(), Zero) {}

template <class Type>
Foam::pulseFixedValueFvPatchField<Type>::pulseFixedValueFvPatchField(
    const fvPatch &p, const DimensionedField<Type, volMesh> &iF,
    const dictionary &dict)
    : fixedValueFvPatchField<Type>(p, iF),
      periodo_(dict.lookup<scalar>("periodo")),
      duracion_(dict.lookup<scalar>("duracion")),
      baseValue_(dict.lookupOrDefault<Type>("baseValue", Zero)),
      pulseValue_("pulseValue", dict, p.size()) {

  fixedValueFvPatchField<Type>::evaluate();

  /*
  // Initialise with the value entry if evaluation is not possible
  fvPatchField<Type>::operator=
  (
      Field<Type>("value", dict, p.size())
  );
  */
}

template <class Type>
Foam::pulseFixedValueFvPatchField<Type>::pulseFixedValueFvPatchField(
    const pulseFixedValueFvPatchField<Type> &ptf, const fvPatch &p,
    const DimensionedField<Type, volMesh> &iF, const fvPatchFieldMapper &mapper)
    : fixedValueFvPatchField<Type>(ptf, p, iF, mapper), periodo_(ptf.periodo_),duracion_(ptf.duracion_),
      baseValue_(ptf.baseValue_), pulseValue_(mapper(ptf.pulseValue_)) {}

template <class Type>
Foam::pulseFixedValueFvPatchField<Type>::pulseFixedValueFvPatchField(
    const pulseFixedValueFvPatchField<Type> &ptf,
    const DimensionedField<Type, volMesh> &iF)
    : fixedValueFvPatchField<Type>(ptf, iF), periodo_(ptf.periodo_),duracion_(ptf.duracion_),
      baseValue_(ptf.baseValue_), pulseValue_(ptf.pulseValue_) {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Type>
void Foam::pulseFixedValueFvPatchField<Type>::autoMap(
    const fvPatchFieldMapper &m) {
  fixedValueFvPatchField<Type>::autoMap(m);
  m(pulseValue_, pulseValue_);
}

      template<class Type>
      void Foam::pulseFixedValueFvPatchField<Type>::rmap
      (
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
      {
  fixedValueFvPatchField<Type>::rmap(ptf, addr);

  const pulseFixedValueFvPatchField<Type> &tiptf =
      refCast<const pulseFixedValueFvPatchField<Type>>(ptf);

  pulseValue_.rmap(tiptf.pulseValue_, addr);
      }


      template<class Type>
      void Foam::pulseFixedValueFvPatchField<Type>::updateCoeffs()
      {
  if (this->updated()) {
    return;
  }

  fixedValueFvPatchField<Type>::operator==(baseValue_ + pulseValue_);

  fixedValueFvPatchField<Type>::updateCoeffs();
      }


      template<class Type>
      void Foam::pulseFixedValueFvPatchField<Type>::write
      (
    Ostream& os
) const
      {
  fvPatchField<Type>::write(os);
  writeEntry(os, "periodo", periodo_);
  writeEntry(os, "duracion", duracion_);
  writeEntry(os, "baseValue", baseValue_);
  writeEntry(os, "pulseValue", pulseValue_);
  writeEntry(os, "value", *this);
      }



// ************************************************************************* //
