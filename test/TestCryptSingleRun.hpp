#ifndef TESTCRYPTSINGLERUN_HPP_
#define TESTCRYPTSINGLERUN_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"

#include "OffLatticeSimulation.hpp"
#include "VolumeTrackingModifier.hpp"

#include "CryptCellsGenerator.hpp"
#include "CellsGenerator.hpp"

#include "CryptCellCycleModel.hpp"

#include "PanethCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellParentId.hpp"

#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellParentIdWriter.hpp"
#include "CellLabelWriter.hpp"
#include "SeparatedCellLabelWriter.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentialAdhesionSpringForce.hpp"
#include "RepulsionForce.hpp"
#include "CellRetainerForce.hpp"
#include "CellLabel.hpp"

#include "CryptGeometryBoundaryCondition3d.hpp"

#include "VariableSeparationCentreBasedDivisionRule.hpp"
#include "NodeBasedCellPopulation.hpp"

#include "CellLocationTrackingModifier.hpp"

#include "SloughingCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"

#include "PetscSetupAndFinalize.hpp"

#include "SmartPointers.hpp"
#include "CommandLineArguments.hpp"

class TestCryptSingleRun : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void Test3DCrypt() throw (Exception)
    {
        double proportion_labeled_cells = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-labeled_ratio");
        double separation_multiplier = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-separation_multiplier");
        
        // double proportion_labeled_cells = 0.5;
        // double separation_multiplier = 3.0;
        double default_separation = 1.0*7.0; // 1 CD


        double end_time = 11000;
        // Crypt Setup
        double cell_radius = 3.5;
        double crypt_length = 70; //70
        double crypt_radius = 8.0/M_PI*6.0;                 // Choosing same dimensions as for halted migration paper
        // For this size domain there are about 75 cells in the bottom hemisphere so this makes about 20% Paneth cells

        unsigned num_paneth_cells = 15;//15;
        unsigned num_stem_cells = 60; //60;
        unsigned num_cells = num_paneth_cells + num_stem_cells;

        double stem_retainer_force_magnitude = 7.5*10;
        double paneth_retainer_force_magnitude = 7.5*10;

        // Cell Model setup
        // 1 - Pedigree
        // 2 - Spatial Wnt
        // 3 - Spatial Wnt at birth
        // 4 - Mutant
        unsigned cell_proliferation_model = 3;  // Optimal from MBOC paper 
        bool wnt_dependent_ccd = true;   // Optimal from MBOC paper
        double param = 0.6; // Optimal fr0m MBOC paper
        double CIparam = 0.9; // Optimal from MBOC paper

        // Create some starter nodes
        std::vector<Node<3>*> nodes;
        for(unsigned node_index= 0;  node_index<num_cells; node_index++)
        {
            double x = crypt_radius/2.0 * sin(node_index*2.0*M_PI/num_cells);
            double y = crypt_radius/2.0 * cos(node_index*2.0*M_PI/num_cells);
            double z = 0.0;
            nodes.push_back(new Node<3>(node_index, false, x, y, z));
        }

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,cell_radius*3.0);

        MAKE_PTR(PanethCellProliferativeType, p_paneth_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        MAKE_PTR(CellLabel,p_label);

        // Create cells
        std::vector<CellPtr> cells;

        CellsGenerator<CryptCellCycleModel, 3> cells_generator;
        
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(),p_transit_type);

        //Change properties of the ccm
        for (unsigned cell_index= 0;  cell_index<cells.size(); cell_index++)
        {

            cells[cell_index]->GetCellData()->SetItem("Radius", cell_radius);
            
            // set up the parent IDs
            MAKE_PTR_ARGS(CellParentId, p_cell_parent_id,(cell_index));
            cells[cell_index]->AddCellProperty(p_cell_parent_id);

            cells[cell_index]->GetCellData()->SetItem("Radius", cell_radius);

            // Initialise the division location on initial cells
            cells[cell_index]->GetCellData()->SetItem("division_location_x",0);
            cells[cell_index]->GetCellData()->SetItem("division_location_y",0);
            cells[cell_index]->GetCellData()->SetItem("division_location_z",0);
            cells[cell_index]->GetCellData()->SetItem("division_time",0);
            cells[cell_index]->GetCellData()->SetItem("parent_id",-1);


            // // Specify CCM
            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetCellProliferationModel(cell_proliferation_model);
            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetIsWntDependentCellCycleDuration(wnt_dependent_ccd);

            // Set some default CCD parameters So total CCM is U[10,14] and (U[22,26] at base if variable)
            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetMDuration(4.0);
            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetSDuration(4.0);
            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetG2Duration(2.0);
            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetTransitCellG1Duration(2.0);  // so total CCM is U[10,14] at threshold
            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetStemCellG1Duration(14.0);  // so total CCM is U[10,14] at base

            //Threshold and Generation specific parameters
            if (cell_proliferation_model ==1 ) // i.e Pedigree dependent
            {
                // Not used as only using Model 3
                NEVER_REACHED;
            }
            else
            {
                dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetWntThreshold(param);
                dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetMaxTransitGenerations(UINT_MAX);
            }

            // Contact Inhibition specific parameters
            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetEquilibriumVolume(M_PI*4.0/3.0*cell_radius*cell_radius*cell_radius);
            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetQuiescentVolumeFraction(CIparam);

            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetLabelledProbability(proportion_labeled_cells);
        }   

        // Make some cells paneth cells
        for(unsigned cell_index= 0;  cell_index<num_paneth_cells; cell_index++)
        {
            unsigned index = cell_index * num_cells/num_paneth_cells;
            if (index > num_cells)
            {
                index = num_cells;
            }
            cells[index]->SetCellProliferativeType(p_paneth_type);
        }

        // Create cell population
        NodeBasedCellPopulation<3> crypt(mesh, cells);
        crypt.SetUseVariableRadii(true);

        // Output data
        crypt.AddCellWriter<CellAgesWriter>();
        crypt.AddCellWriter<CellVolumesWriter>();
        crypt.AddCellWriter<CellProliferativeTypesWriter>();
        crypt.AddCellWriter<CellMutationStatesWriter>();
        crypt.AddCellWriter<CellIdWriter>();
        crypt.AddCellWriter<CellParentIdWriter>();
        crypt.AddCellWriter<CellLabelWriter>();
        crypt.AddCellWriter<SeparatedCellLabelWriter>();
        
        crypt.AddPopulationWriter<NodeVelocityWriter>();
        crypt.SetAbsoluteMovementThreshold(50.0);

        // Set the division rule for our population to be the random direction division rule
        boost::shared_ptr<VariableSeparationCentreBasedDivisionRule<3,3> > p_division_rule(new VariableSeparationCentreBasedDivisionRule<3,3>());
        p_division_rule->SetCellSeparation(default_separation);
        p_division_rule->SetLabeledSeparationMultiplier(separation_multiplier);
        crypt.SetCentreBasedDivisionRule(p_division_rule);

        // Create an instance of a Wnt concentration NOTE DO THIS BEFORE THE SIMULATION OTHERWISE CELLS CANT INITIALISE
        WntConcentration<3>::Instance()->SetType(LINEAR);
        WntConcentration<3>::Instance()->SetCellPopulation(crypt);
        WntConcentration<3>::Instance()->SetCryptLength(crypt_length);

        // Create a contact inhibition simulator
        OffLatticeSimulation<3> simulator(crypt);

        simulator.SetOutputDivisionLocations(true);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(end_time);

        // Add Volume Tracking Modifier
        MAKE_PTR(VolumeTrackingModifier<3>, p_volume_modifier);
        simulator.AddSimulationModifier(p_volume_modifier);

        // Add Modoifier to allow cells to track locations
        MAKE_PTR(CellLocationTrackingModifier<3>,p_location_modifier);
        simulator.AddSimulationModifier(p_location_modifier);

        //Create output directory
        std::stringstream out;
        out << "Ratio_"<< proportion_labeled_cells << "/Separation_" << separation_multiplier;
        std::string output_directory = "3dCryptCellSeparation" +  out.str();
        simulator.SetOutputDirectory(output_directory);

        // Create a force law and pass it to the simulation
        MAKE_PTR(DifferentialAdhesionSpringForce<3>, p_force);
        p_force->SetMeinekeSpringStiffness(30.0); //normally 15.0 but 30 in all CellBased Papers;
        p_force->SetCutOffLength(cell_radius*3.0); // (1.5 * celldiameter)
        simulator.AddForce(p_force);
                
        // Apply a retainer to keep stem and paneth cells at the base of the crypt
        MAKE_PTR(CellRetainerForce<3>, p_retainer_force);
        p_retainer_force->SetStemCellForceMagnitudeParameter(stem_retainer_force_magnitude);
        p_retainer_force->SetPanethCellForceMagnitudeParameter(paneth_retainer_force_magnitude);
        simulator.AddForce(p_retainer_force);

        //Apply a boundary condition to represent a 3d crypt
        MAKE_PTR_ARGS(CryptGeometryBoundaryCondition3d<3>, p_boundary_condition, (&crypt, 0.0));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer,(&crypt, crypt_length*unit_vector<double>(3,2), unit_vector<double>(3,2)));
        simulator.AddCellKiller(p_cell_killer);

        // Run simulation
        simulator.Solve();
    }
};

#endif /*TESTCRYPTSINGLERUN_HPP_*/
