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

#include "HoneycombMeshGenerator.hpp"
#include "Cylindrical2dNodesOnlyMesh.hpp"

#include "CryptCellCycleModel.hpp"
#include "NewSimpleWntCellCycleModel.hpp"

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

#include "PlaneBasedCellKiller.hpp"
#include "PlaneBoundaryCondition.hpp"

#include "PetscSetupAndFinalize.hpp"

#include "SmartPointers.hpp"
#include "CommandLineArguments.hpp"

#include "Debug.hpp"

class Test2dCryptSingleRun : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void Test2DCrypt() throw (Exception)
    {
    	double proportion_labeled_cells = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-labeled_ratio");
    	double separation_multiplier = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-separation_multiplier");
    	
    	//double proportion_labeled_cells = 0.5;
		//double separation_multiplier = 3.0;
		double default_separation = 1.0;


    	double end_time = 2200;
        // Crypt Setup
    	double crypt_width = 16; //16
        double crypt_length = 12; // Dimensions from Cell Migration Plos One

    	// Cell Model setup
		unsigned cell_proliferation_model = 2; // 3  // Optimal from MBOC paper Spatial Wnt at birth
		bool wnt_dependent_ccd = false;   // Optimal from MBOC paper is true!!!
        double param = 0.5; // 
        double CIparam = 0.0; // No CI at moment all cells proliferate

        // Create a simple mesh
        HoneycombMeshGenerator generator(crypt_width, crypt_length, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        double cut_off_length = 1.5; //this is the default

        // Convert this to a Cylindrical2dNodesOnlyMesh
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(crypt_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh,2.0); // So factor of crypt_width

          
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);

		MAKE_PTR(CellLabel,p_label);

		// Create cells
		std::vector<CellPtr> cells;
		CellsGenerator<CryptCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(),p_transit_type);

		//Change properties of the ccm
		for (unsigned cell_index= 0;  cell_index<p_mesh->GetNumNodes(); cell_index++)
		{
			// set up the parent IDs
			MAKE_PTR_ARGS(CellParentId, p_cell_parent_id,(cell_index));
			cells[cell_index]->AddCellProperty(p_cell_parent_id);

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


			// SetWnt Specific parameters
		    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetWntThreshold(param);
		    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetMaxTransitGenerations(UINT_MAX);
			    
            // Contact Inhibition specific parameters
            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetEquilibriumVolume(sqrt(3.0)/2.0);
            dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetQuiescentVolumeFraction(CIparam);
		}

		// Create cell population
		NodeBasedCellPopulation<2> crypt(*p_mesh, cells);

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
		boost::shared_ptr<VariableSeparationCentreBasedDivisionRule<2,2> > p_division_rule(new VariableSeparationCentreBasedDivisionRule<2,2>());
		p_division_rule->SetCellSeparation(default_separation);
		p_division_rule->SetLabeledSeparationMultiplier(separation_multiplier);
		crypt.SetCentreBasedDivisionRule(p_division_rule);

		// Create an instance of a Wnt concentratizon NOTE DO THIS BEFORE THE SIMULATION OTHERWISE CELLS CANT INITIALISE
		WntConcentration<2>::Instance()->SetType(LINEAR);
		WntConcentration<2>::Instance()->SetCellPopulation(crypt);
		WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

		// Create a contact inhibition simulator
		OffLatticeSimulation<2> simulator(crypt);

		simulator.SetOutputDivisionLocations(true);
		simulator.SetDt(1.0/200.0);
		simulator.SetSamplingTimestepMultiple(200);
		simulator.SetEndTime(end_time);

		// Add Volume Tracking Modifier
		MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
		simulator.AddSimulationModifier(p_volume_modifier);

		// Add Modoifier to allow cells to track locations
		MAKE_PTR(CellLocationTrackingModifier<2>,p_location_modifier);
		simulator.AddSimulationModifier(p_location_modifier);

		//Create output directory
		std::stringstream out;
		out << "Ratio_"<< proportion_labeled_cells << "/Separation_" << separation_multiplier;
		std::string output_directory = "2dCryptCellSeparation/" +  out.str();
		simulator.SetOutputDirectory(output_directory);

		// Create a force law and pass it to the simulation
		MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
		p_force->SetMeinekeSpringStiffness(50.0); //normally 15.0 but 30 in all CellBased Papers;
		p_force->SetCutOffLength(1.5);
		simulator.AddForce(p_force);
		
        // Solid base boundary condition
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bcs, (&crypt, zero_vector<double>(2), -unit_vector<double>(2,1)));
        p_bcs->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bcs);

		// Create cell killer and pass in to crypt simulation
		MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_cell_killer,(&crypt, crypt_length*unit_vector<double>(2,1), unit_vector<double>(2,1)));
		simulator.AddCellKiller(p_cell_killer);

		// Run simulation
		simulator.Solve();

        // Clear memory
        delete p_mesh;
	}
};

#endif /*TESTCRYPTSINGLERUN_HPP_*/
