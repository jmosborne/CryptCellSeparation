/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef CRYPTGEOMETRYBOUNDARYCONDITION3D_HPP_
#define CRYPTGEOMETRYBOUNDARYCONDITION3D_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A spherical cell population boundary condition class, which restricts nodes to lie
 * on the surface of a deformed test tube in the domain.
 */
template<unsigned DIM>
class CryptGeometryBoundaryCondition3d : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:

    /** The maximum distance from the surface that cells may be. */
    double mMaximumDistance;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
        archive & mMaximumDistance;
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param radius the radius of the base
     * @param distance the maximum distance from the surface that cells may be (defaults to 1e-5)
     */
    CryptGeometryBoundaryCondition3d(AbstractCellPopulation<DIM>* pCellPopulation, double distance=1e-5);

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CryptGeometryBoundaryCondition3d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptGeometryBoundaryCondition3d.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CryptGeometryBoundaryCondition3d<DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
}

/**
 * De-serialize constructor parameters and initialize a CryptGeometryBoundaryCondition3d.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CryptGeometryBoundaryCondition3d<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)CryptGeometryBoundaryCondition3d<DIM>(p_cell_population);
}
}
} // namespace ...

#endif /*CRYPTGEOMETRYBOUNDARYCONDITION3D_HPP_*/
