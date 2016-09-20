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

#include "CryptGeometryBoundaryCondition3d.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"

/* This is a boundary condition which is determined from Paul Appleton's crypt data - we have fit a polynomial
 * to the data to generate a function for the crypt shape. An ellipse is positioned to generate the crypt base.
 */


template<unsigned DIM>
CryptGeometryBoundaryCondition3d<DIM>::CryptGeometryBoundaryCondition3d(AbstractCellPopulation<DIM>* pCellPopulation,
                                                                      double distance)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
      mMaximumDistance(distance)
{
    assert(mMaximumDistance >= 0.0);

    if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation) == NULL)
    {
        EXCEPTION("A NodeBasedCellPopulation must be used with this boundary condition object.");
    }
    if (DIM < 3)
    {
        EXCEPTION("This boundary condition is only implemented in 3D.");
    }
}

template<unsigned DIM>
void CryptGeometryBoundaryCondition3d<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
	// Set the radius of the crypt base appropriately (an ellipse)
	double major_radius = 16.6968;
	double minor_radius = 15.3973;

	double factor = 1.0;	// Using this to get it roughly to a similar size as the previous crypt shape

	major_radius /= factor;
	minor_radius /= factor;

	// polynomial coefficients: p = p_1x^5 + p_2x^4 + p_3x^3 + p_4x^2 + p_5x + p_6
	double p_1 = 8.365027421137118e-07;
	double p_2 = -1.612884286709950e-04;
	double p_3 = 1.169561186127581e-02;
	double p_4 = -3.912727534949710e-01;
	double p_5 = 5.850759485536169e+00;
	double p_6 = -1.497897551314556e+01
;

	// Iterate over the cell population
	    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
	         cell_iter != this->mpCellPopulation->End();
	         ++cell_iter)
	    {
	        c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

	        c_vector<double,DIM> ellipse_centre = zero_vector<double>(DIM);
	        ellipse_centre[DIM-1] = minor_radius;	// the shorter of the two radii in the ellipse

	        // Where the cell originally sat along the z axis (will be used to put it back to its original height later)
	        double original_height = cell_location[DIM-1];

	        // The target radius depends on each node, and is defined separately for those cells in the crypt base,
	        // and along the crypt walls
	        double target_radius = 0.0;

	        // Those cells that sit at the crypt base
	        if (cell_location[DIM-1] < minor_radius)
	        {
	        	double target_x_coordinate = sqrt((pow(major_radius,2)*pow(minor_radius,2)*pow(cell_location[DIM-3],2)) /
	        			(pow(minor_radius,2)*pow(cell_location[DIM-3],2) + pow(major_radius,2)*pow(cell_location[DIM-2],2)));
	        	double target_y_coordinate = sqrt((pow(major_radius,2)*pow(minor_radius,2)*pow(cell_location[DIM-2],2)) /
	        			(pow(minor_radius,2)*pow(cell_location[DIM-3],2) + pow(major_radius,2)*pow(cell_location[DIM-2],2)));

	        	target_radius = sqrt(pow(target_x_coordinate,2) + pow(target_y_coordinate,2));
	        }

	        // For those cells which sit above the elliptic base of the crypt
	        if (cell_location[DIM-1] >= minor_radius)
	        {
	            double z = cell_location[DIM-1]*factor; // To avoid the base of the function where it goes negative

	            // Get the target radius for this cell (which sits higher in the crypt)
	        	target_radius = p_1*z*z*z*z*z + p_2*z*z*z*z + p_3*z*z*z + p_4*z*z + p_5*z + p_6; // interpolated crypt shape
	        	//target_radius /= factor;
//	        	PRINT_3_VARIABLES(p_1, p_2, p_3);
//
//	        	PRINT_3_VARIABLES(p_4, p_5, p_6);
//	        	PRINT_2_VARIABLES(z*z, z);
//
//	        	PRINT_3_VARIABLES(p_1*z*z*z*z*z, p_2*z*z*z*z, p_3*z*z*z);
//	        	PRINT_3_VARIABLES(p_4*z*z, p_5*z, p_6);
//	            PRINT_2_VARIABLES(z, target_radius);
//TRACE("NewLine")
	            cell_location[DIM-1] = minor_radius;
	        }

	        // Find the radial distance between this cell and the surface
	        double radius = norm_2(cell_location - ellipse_centre);

	        // If the cell is too far from the surface of the ellipsoid / upper crypt
	        if (fabs(radius - target_radius) > mMaximumDistance)
	        {
	            // ...move the cell back onto the surface

	            c_vector<double, DIM> location_on_surface = target_radius*(cell_location - ellipse_centre)/radius + ellipse_centre;

	            if (original_height >= minor_radius)
	            {
	                assert(fabs(location_on_surface[DIM-1]-minor_radius)<1e-8);
	                location_on_surface[DIM-1] = original_height;
	            }

	            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
	            Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

	            if ( location_on_surface[DIM-1] < 0 )
	            {
	                //PRINT_2_VARIABLES(node_index,location_on_surface[DIM-1]);
	                location_on_surface[DIM-1] = 0.0;
	            }

	            p_node->rGetModifiableLocation() = location_on_surface;

	        }
	    }
}

template<unsigned DIM>
bool CryptGeometryBoundaryCondition3d<DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    // Iterate over the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
//        c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
//
//        c_vector<double,DIM> sphere_centre = zero_vector<double>(DIM);
//        sphere_centre[DIM-1] = mRadiusOfBase;
//        //double target_radius = mRadiusOfBase;
//        //double original_height = cell_location[DIM-1];
//
//        if (cell_location[DIM-1] > mRadiusOfBase)
//        {
//            //double z = cell_location[DIM-1] - mRadiusOfBase;
//            //target_radius *= (1.0525 - 0.05*(tanh(0.5*(z-5.0)) - 2*tanh(1*(z-9.0))));
//            cell_location[DIM-1]=mRadiusOfBase;
//        }
//
//        // Find the radial distance between this cell and the surface
//        double radius = norm_2(cell_location - sphere_centre);
//
//        // If the cell is too far from the surface of the sphere...
//        if (fabs(radius - mRadiusOfBase) > mMaximumDistance+1e-12)
//        {
//            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
//            condition_satisfied = false;
//            PRINT_4_VARIABLES(node_index, cell_location[0],cell_location[1],cell_location[2]);
//            PRINT_2_VARIABLES(radius,mRadiusOfBase);
//        }
    }


    return condition_satisfied;
}

template<unsigned DIM>
void CryptGeometryBoundaryCondition3d<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MaximumDistance>" << mMaximumDistance << "</MaximumDistance>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CryptGeometryBoundaryCondition3d<1>;
template class CryptGeometryBoundaryCondition3d<2>;
template class CryptGeometryBoundaryCondition3d<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CryptGeometryBoundaryCondition3d)
