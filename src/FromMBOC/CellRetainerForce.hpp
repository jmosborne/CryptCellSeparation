/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef CELLRETAINERFORCE_HPP_
#define CELLRETAINERFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"

/**
 * A force class to retain stem and paneth cells in the base of the crypt.
 */
template<unsigned DIM>
class CellRetainerForce : public AbstractForce<DIM>
{
private :

    /**
     * Force applied to stem cells
     */
    double mStemCellForceMagnitudeParameter;

    /**
     * Force applied to paneth cells
     */
    double mPanethCellForceMagnitudeParameter;

    /**
     * Archiving.
     */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mStemCellForceMagnitudeParameter;
        archive & mPanethCellForceMagnitudeParameter;
    }

public :

    /**
     * Constructor.
     */
    CellRetainerForce();

    /**
     * Destructor.
     */
    ~CellRetainerForce();

    /**
      * Set the magnitude of the retainer force.
      *
      * @param stemCellForceMagnitudeParameter the new parameter
      */
    void SetStemCellForceMagnitudeParameter(double stemCellForceMagnitudeParameter);

    /**
     * Get the StemCellForceMagnitudeParameter.
     *
     * @return mStemCellForceMagnitudeParameter
     */
    double GetStemCellForceMagnitudeParameter();

    /**
      * Set the magnitude of the retainer force.
      *
      * @param panetCellForceMagnitudeParameter the new parameter
      */
    void SetPanethCellForceMagnitudeParameter(double panethCellForceMagnitudeParameter);

    /**
     * Get the PanethCellForceMagnitudeParameter.
     *
     * @return mPanethCellForceMagnitudeParameter
     */
    double GetPanethCellForceMagnitudeParameter();

     /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the tissue
     *
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellRetainerForce)

#endif /*CELLRETAINERFORCE_HPP_*/
