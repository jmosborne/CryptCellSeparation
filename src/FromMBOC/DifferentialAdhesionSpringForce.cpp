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

#include "DifferentialAdhesionSpringForce.hpp"
#include "PanethCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "Debug.hpp"

template<unsigned DIM>
DifferentialAdhesionSpringForce<DIM>::DifferentialAdhesionSpringForce()
   : GeneralisedLinearSpringForce<DIM>()
{
}

template<unsigned DIM>
double DifferentialAdhesionSpringForce<DIM>::VariableSpringConstantMultiplicationFactor(
    unsigned nodeAGlobalIndex,
    unsigned nodeBGlobalIndex,
    AbstractCellPopulation<DIM>& rCellPopulation,
    bool isCloserThanRestLength)
{

    double multiplication_factor = 1.0;

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    if ((p_cell_A->GetCellProliferativeType()->template IsType<PanethCellProliferativeType>())
        && (p_cell_B->GetCellProliferativeType()->template IsType<PanethCellProliferativeType>())
        && !isCloserThanRestLength)
    {
        multiplication_factor =  0.0;
    }


    return multiplication_factor;
}


// The only diference here is the alpha parameter in the node based attraction!!!! (AND SPACE/ELEMEMENT_DIM => DIM and some this's)
template<unsigned DIM>
c_vector<double, DIM> DifferentialAdhesionSpringForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, DIM> node_a_location = p_node_a->rGetLocation();
    c_vector<double, DIM> node_b_location = p_node_b->rGetLocation();

    // Get the node radii for a NodeBasedCellPopulation
    double node_a_radius=0.0;
    double node_b_radius=0.0;

    if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        node_a_radius = p_node_a->GetRadius();
        node_b_radius = p_node_b->GetRadius();
    }

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, DIM> unit_difference;
    /*
     * We use the mesh method GetVectorFromAtoB() to compute the direction of the
     * unit vector along the line joining the two nodes, rather than simply subtract
     * their positions, because this method can be overloaded (e.g. to enforce a
     * periodic boundary in Cylindrical2dMesh).
     */
    unit_difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);


    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    /*
     * If mUseCutOffLength has been set, then there is zero force between
     * two nodes located a distance apart greater than mMechanicsCutOffLength in AbstractTwoBodyInteractionForce.
     */
    if (this->mUseCutOffLength)
    {
        if (distance_between_nodes >= this->GetCutOffLength())
        {
            return zero_vector<double>(DIM); // c_vector<double,DIM>() is not guaranteed to be fresh memory
        }
    }

    /*
     * Calculate the rest length of the spring connecting the two nodes with a default
     * value of 1.0.
     */
    double rest_length_final = 1.0;

    if (dynamic_cast<MeshBasedCellPopulation<DIM,DIM>*>(&rCellPopulation))
    {
        rest_length_final = static_cast<MeshBasedCellPopulation<DIM,DIM>*>(&rCellPopulation)->GetRestLength(nodeAGlobalIndex, nodeBGlobalIndex);
    }
    else if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        assert(node_a_radius > 0 && node_b_radius > 0);
        rest_length_final = node_a_radius+node_b_radius;
    }

    double rest_length = rest_length_final;

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    double ageA = p_cell_A->GetAge();
    double ageB = p_cell_B->GetAge();

    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));

    /*
     * If the cells are both newly divided, then the rest length of the spring
     * connecting them grows linearly with time, until 1 hour after division.
     */
    if (ageA < this->mMeinekeSpringGrowthDuration && ageB < this->mMeinekeSpringGrowthDuration)
    {
        AbstractCentreBasedCellPopulation<DIM,DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<DIM,DIM>*>(&rCellPopulation);

        std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

        if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
        {
            // Spring rest length increases from a small value to the normal rest length over 1 hour
            double lambda = this->mMeinekeDivisionRestingSpringLength;
            rest_length = lambda + (rest_length_final - lambda) * ageA/this->mMeinekeSpringGrowthDuration;
        }
        if (ageA + SimulationTime::Instance()->GetTimeStep() >= this->mMeinekeSpringGrowthDuration)
        {
            // This spring is about to go out of scope
            p_static_cast_cell_population->UnmarkSpring(cell_pair);
        }
    }

    /*
     * For apoptosis, progressively reduce the radius of the cell
     */
    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        assert(node_a_radius > 0 && node_b_radius > 0);
        a_rest_length = (node_a_radius/(node_a_radius+node_b_radius))*rest_length;
        b_rest_length = (node_b_radius/(node_a_radius+node_b_radius))*rest_length;
    }

    /*
     * If either of the cells has begun apoptosis, then the length of the spring
     * connecting them decreases linearly with time.
     */
    if (p_cell_A->HasApoptosisBegun())
    {
        double time_until_death_a = p_cell_A->GetTimeUntilDeath();
        a_rest_length = a_rest_length * time_until_death_a / p_cell_A->GetApoptosisTime();
    }
    if (p_cell_B->HasApoptosisBegun())
    {
        double time_until_death_b = p_cell_B->GetTimeUntilDeath();
        b_rest_length = b_rest_length * time_until_death_b / p_cell_B->GetApoptosisTime();
    }

    rest_length = a_rest_length + b_rest_length;
    //assert(rest_length <= 1.0+1e-12); ///\todo #1884 Magic number: would "<= 1.0" do?


    // Although in this class the 'spring constant' is a constant parameter, in
    // subclasses it can depend on properties of each of the cells
    double overlap = distance_between_nodes - rest_length;
    bool is_closer_than_rest_length = (overlap <= 0);
    double multiplication_factor = VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex, nodeBGlobalIndex, rCellPopulation, is_closer_than_rest_length);
    double spring_stiffness = this->mMeinekeSpringStiffness;

    if (dynamic_cast<MeshBasedCellPopulation<DIM,DIM>*>(&rCellPopulation))
    {
        return multiplication_factor * spring_stiffness * unit_difference * overlap;
    }
    else
    {
        // A reasonably stable simple force law
        if (is_closer_than_rest_length) //overlap is negative
        {
            //log(x+1) is undefined for x<=-1
            assert(overlap > -rest_length_final);
            c_vector<double, DIM> temp = multiplication_factor*spring_stiffness * unit_difference * rest_length_final* log(1.0 + overlap/rest_length_final);
            return temp;
        }
        else
        {
            double alpha = 4;
            c_vector<double, DIM> temp = multiplication_factor*spring_stiffness * unit_difference * overlap * exp(-alpha * overlap/rest_length_final);
            return temp;
        }
    }
}

template<unsigned DIM>
void DifferentialAdhesionSpringForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call direct parent class
    GeneralisedLinearSpringForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DifferentialAdhesionSpringForce<1>;
template class DifferentialAdhesionSpringForce<2>;
template class DifferentialAdhesionSpringForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DifferentialAdhesionSpringForce)
