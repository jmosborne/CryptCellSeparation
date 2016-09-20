/*

Copyright (c) 2005-2013, University of Oxford.
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

#include "CryptCellCycleModel.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "PanethCellProliferativeType.hpp"
#include "WntConcentration.hpp"
#include "CellLabel.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SeparatedCellLabel.hpp"
#include "CellParentId.hpp"
#include "Debug.hpp"


CryptCellCycleModel::CryptCellCycleModel()
: AbstractSimpleGenerationalCellCycleModel(),
  mCellProliferationModel(UNSIGNED_UNSET),
  mIsWntDependentCellCycleDuration(false),
  mWntThreshold(DOUBLE_UNSET),
  mInitialWntLevel(1.0), // Initially all cells will be stem like cells
  mQuiescentVolumeFraction(DOUBLE_UNSET),
  mEquilibriumVolume(DOUBLE_UNSET),
  mCurrentQuiescentOnsetTime(SimulationTime::Instance()->GetTime()),
  mCurrentQuiescentDuration(0.0),
  mLabelledProbability(0.5)
{
}

CryptCellCycleModel::CryptCellCycleModel(const CryptCellCycleModel& rModel)
: AbstractSimpleGenerationalCellCycleModel(rModel),
  mCellProliferationModel(rModel.mCellProliferationModel),
  mIsWntDependentCellCycleDuration(rModel.mIsWntDependentCellCycleDuration),
  mWntThreshold(rModel.mWntThreshold),
  mInitialWntLevel(rModel.mInitialWntLevel), 
  mQuiescentVolumeFraction(rModel.mQuiescentVolumeFraction),
  mEquilibriumVolume(rModel.mEquilibriumVolume),
  mCurrentQuiescentOnsetTime(rModel.mCurrentQuiescentOnsetTime),
  mCurrentQuiescentDuration(rModel.mCurrentQuiescentDuration),
  mLabelledProbability(rModel.mLabelledProbability)
{
}

void CryptCellCycleModel::UpdateCellCyclePhase()
{

    // Paneth Cells remain differentiated under all models
    if (mpCell->GetCellProliferativeType()->IsType<PanethCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else // Everything but paneth
    {
        // Implement contact inhibition

        // Get cell volume
        double cell_volume = mpCell->GetCellData()->GetItem("volume");

        // Removes the cell label
        mpCell->RemoveCellProperty<CellLabel>();

        if (mCurrentCellCyclePhase == G_ONE_PHASE)
        {
               // Update G1 duration based on cell volume
               double dt = SimulationTime::Instance()->GetTimeStep();


               double quiescent_volume = mEquilibriumVolume * mQuiescentVolumeFraction;

               if (cell_volume < quiescent_volume)
               {
                       // Update the duration of the current period of contact inhibition.
                       mCurrentQuiescentDuration = SimulationTime::Instance()->GetTime() - mCurrentQuiescentOnsetTime;
                       mG1Duration += dt;
                       mpCell->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellLabel>());
               }
                else
               {
                       // Reset the cell's quiescent duration and update the time at which the onset of quiescent occurs
                       mCurrentQuiescentDuration = 0.0;
                       mCurrentQuiescentOnsetTime = SimulationTime::Instance()->GetTime();
               }
        }
      

        if(mCellProliferationModel == 1) // Pedigree model
        {
            AbstractSimpleGenerationalCellCycleModel::UpdateCellCyclePhase();
        }
        else if(mCellProliferationModel == 2 || mCellProliferationModel == 3) // Spatial Wnt and Spatial Wnt on Birth
        {
            // The cell can divide if the Wnt concentration >= GetWntThreshold()

            double wnt_level;

            // Only difference between 2 models is where they take wnt level from
            if(mCellProliferationModel == 2)
            {
                wnt_level= GetWntLevel();
            }
            else
            {
                assert(mCellProliferationModel == 3);
                wnt_level = mInitialWntLevel;
             }

            // Set the cell type to TransitCellProliferativeType if the Wnt stimulus exceeds wnt_division_threshold if not set it to Differentiated Type
            if (wnt_level >= GetWntThreshold())
            {
                boost::shared_ptr<AbstractCellProperty> p_transit_type =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
                mpCell->SetCellProliferativeType(p_transit_type);
            }
            else
            {
                // The cell is set to have DifferentiatedCellProliferativeType and so in G0 phase
                boost::shared_ptr<AbstractCellProperty> p_diff_type =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
                mpCell->SetCellProliferativeType(p_diff_type);
            }
            AbstractSimplePhaseBasedCellCycleModel::UpdateCellCyclePhase();
        }
        else if(mCellProliferationModel == 4) // mutant
		{
        	NEVER_REACHED;
        	// All Cells are proliferative
        	boost::shared_ptr<AbstractCellProperty> p_transit_type =
        	    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        	mpCell->SetCellProliferativeType(p_transit_type);

        	AbstractSimplePhaseBasedCellCycleModel::UpdateCellCyclePhase();
		}
        else
        {
            NEVER_REACHED;
        }
    }
	//mpCell->GetCellData()->SetItem("G1Duration", mG1Duration);
}

AbstractCellCycleModel* CryptCellCycleModel::CreateCellCycleModel()
{
    return new CryptCellCycleModel(*this);
}

void CryptCellCycleModel::SetG1Duration()
{
//    TRACE("SetG1Duration");
//    PRINT_VARIABLE(SimulationTime::Instance()->GetTime());

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    assert(mpCell != NULL);

    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() ||
        mpCell->GetCellProliferativeType()->IsType<PanethCellProliferativeType>() )
    {
        mG1Duration = DBL_MAX;
    }
    else if(mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        // Only pedigree model includes stem cells
        //assert(mCellProliferationModel == 1); Can't do this as at the first step there are stem cells.

        double mean_g1_duration;

        if(mIsWntDependentCellCycleDuration)
        {
            mean_g1_duration = GetStemCellG1Duration(); //so U[22,26] for default parameters
        }
        else
        {
            mean_g1_duration = GetTransitCellG1Duration(); //so U[10,14] for default parameters
        }
        mG1Duration = mean_g1_duration  -2.0  + 4*p_gen->ranf();
    }
    else if(mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        double mean_g1_duration;
        if(mCellProliferationModel == 1)
        {
            mean_g1_duration = GetTransitCellG1Duration(); //so U[10,14] for default parameters
        }
        else if (mCellProliferationModel == 4)
        {
        	NEVER_REACHED;
        	mean_g1_duration = GetTransitCellG1Duration(); //so U[10,14] for default parameters
        	//Note mutant is independent of WNT Dependent CCD
        }
        else
        {
            // Should be one of the spatial models both have the same ccd
            assert(mCellProliferationModel == 2 || mCellProliferationModel == 3);

            if(mIsWntDependentCellCycleDuration)
            {
                double mean_g1_duration_at_base = GetStemCellG1Duration(); //so U[22,26] for default parameters
                double mean_g1_duration_at_wnt_threshold = GetTransitCellG1Duration(); //so U[10,14] for default parameters

                double wnt_level = mInitialWntLevel; // Note GetWntLevel() wont work as at this point the cell has no spatial location so we use this cached version;

                if (mWntThreshold == 1.0)
                {
                    mean_g1_duration = GetTransitCellG1Duration();
                }
                else
                {
                    mean_g1_duration = (mean_g1_duration_at_base - mean_g1_duration_at_wnt_threshold)/(1.0 - GetWntThreshold()) * (wnt_level-1.0) + mean_g1_duration_at_base;   //so between U[10,14] and U[22,26] for default parameters
                    //PRINT_4_VARIABLES(mean_g1_duration_at_base,mean_g1_duration_at_wnt_threshold, GetWntThreshold(),wnt_level);
                    //PRINT_VARIABLE(mean_g1_duration);
                }
            }
            else
            {
                mean_g1_duration = GetTransitCellG1Duration(); //so U[10,14] for default parameters
            }
        }

        mG1Duration = mean_g1_duration -2.0 + 4*p_gen->ranf();
    }
    else
    {
        NEVER_REACHED;
    }
    //PRINT_3_VARIABLES(mG1Duration,mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>(),mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>());

    //mpCell->GetCellData()->SetItem("InitialG1Duration", mG1Duration);
}


double CryptCellCycleModel::GetWntLevel()
{
    assert(mpCell != NULL);
    double level = 0;

    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        default:
            NEVER_REACHED;
    }

    return level;
}


void CryptCellCycleModel::ResetForDivision()
{

    // This is the new code to track parent IDs and have labeled and non labeled Cells //
    // Reset the cell parent ID
    CellPropertyCollection cell_parent_id_collection = mpCell->rGetCellPropertyCollection().GetPropertiesType<CellParentId>();
    assert(cell_parent_id_collection.GetSize() == 1);
    boost::shared_ptr<CellParentId> p_cell_parent_id = boost::static_pointer_cast<CellParentId>(cell_parent_id_collection.GetProperty());
    p_cell_parent_id->SetParentId(mpCell->GetCellId());

    switch (mDimension)
    {
        case 3:
        {
            double cell_location_x = mpCell->GetCellData()->GetItem("cell_location_x");
            double cell_location_y = mpCell->GetCellData()->GetItem("cell_location_y");
            double cell_location_z = mpCell->GetCellData()->GetItem("cell_location_z");

            mpCell->GetCellData()->SetItem("division_location_x",cell_location_x);
            mpCell->GetCellData()->SetItem("division_location_y",cell_location_y);
            mpCell->GetCellData()->SetItem("division_location_z",cell_location_z);
            break;
        }        
        default:
        NEVER_REACHED;
    }

    mpCell->GetCellData()->SetItem("division_time",SimulationTime::Instance()->GetTime());
   
    // Now reset whether cell is labeled or not.

    // Removes the cell label
    mpCell->RemoveCellProperty<SeparatedCellLabel>();


    if (RandomNumberGenerator::Instance()->ranf() < mLabelledProbability)
    {
      //TRACE("HELLO");
        boost::shared_ptr<AbstractCellProperty> p_label = mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<SeparatedCellLabel>();
        mpCell->AddCellProperty(p_label);
    }
    //////////////////////////////////////////////////////////////////////////////


    if(mCellProliferationModel==1)// Pedegrie Model
    {
        AbstractSimpleGenerationalCellCycleModel::ResetForDivision();
    }
    else // Spatial Model
    {
    	if (mCellProliferationModel==3)// Wnt on birth Model
    	{
    	    // Reset the initial Wnt level in the current ccm
    	    mInitialWntLevel = GetWntLevel();
    	    // Its copied it in the new ccm on creation
    	}

        AbstractSimplePhaseBasedCellCycleModel::ResetForDivision();
    }
}


void CryptCellCycleModel::InitialiseDaughterCell()
{



	// Call appropriate parent class
	if(mCellProliferationModel == 1) // Pedigree model
    {
        AbstractSimpleGenerationalCellCycleModel::InitialiseDaughterCell();
    }
	else 
	{
		AbstractSimplePhaseBasedCellCycleModel::InitialiseDaughterCell();
	}
}


/*
 * Get Set Methods for CCM options
 */
void CryptCellCycleModel::SetCellProliferationModel(unsigned cellProliferationModel)
{
    mCellProliferationModel = cellProliferationModel;
}

unsigned CryptCellCycleModel::GetCellProliferationModel()
{
    return mCellProliferationModel;
}

void CryptCellCycleModel::SetIsWntDependentCellCycleDuration(bool isWntDependentCellCycleDuration)
{
    mIsWntDependentCellCycleDuration = isWntDependentCellCycleDuration;
}

bool CryptCellCycleModel::GetIsWntDependentCellCycleDuration()
{
    return mIsWntDependentCellCycleDuration;
}

double CryptCellCycleModel::GetWntThreshold()
{
		return mWntThreshold;
}

void CryptCellCycleModel::SetWntThreshold(double wntTransitThreshold)
{
    assert(wntTransitThreshold <= 1.0);
    assert(wntTransitThreshold >= 0.0);
    mWntThreshold = wntTransitThreshold;
}

void CryptCellCycleModel::SetInitialWntLevel(double initialWntLevel)
{
    mInitialWntLevel = initialWntLevel;
}

void CryptCellCycleModel::SetQuiescentVolumeFraction(double quiescentVolumeFraction)
{
    mQuiescentVolumeFraction = quiescentVolumeFraction;
}

double CryptCellCycleModel::GetQuiescentVolumeFraction()
{
    return mQuiescentVolumeFraction;
}

void CryptCellCycleModel::SetEquilibriumVolume(double equilibriumVolume)
{
    mEquilibriumVolume = equilibriumVolume;
}

double CryptCellCycleModel::GetEquilibriumVolume()
{
    return mEquilibriumVolume;
}

void CryptCellCycleModel::SetCurrentQuiescentDuration(double currentQuiescentDuration)
{
    mCurrentQuiescentDuration = currentQuiescentDuration;
}

double CryptCellCycleModel::GetCurrentQuiescentDuration()
{
    return mCurrentQuiescentDuration;
}

void CryptCellCycleModel::SetCurrentQuiescentOnsetTime(double currentQuiescentOnsetTime)
{
    mCurrentQuiescentOnsetTime = currentQuiescentOnsetTime;
}

double CryptCellCycleModel::GetCurrentQuiescentOnsetTime()
{
    return mCurrentQuiescentOnsetTime;
}

double CryptCellCycleModel::GetLabelledProbability() const
{
    return mLabelledProbability;
}

void CryptCellCycleModel::SetLabelledProbability(double labelledProbability)
{
    assert(labelledProbability <= 1.0);
    assert(labelledProbability >= 0.0);
    mLabelledProbability = labelledProbability;
}
void CryptCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellProliferationModel>" << mCellProliferationModel << "</CellProliferationModel>\n";
    *rParamsFile << "\t\t\t<IsWntDependentCellCycleDuration>" << mIsWntDependentCellCycleDuration << "</IsWntDependentCellCycleDuration>\n";

    *rParamsFile << "\t\t\t<QuiescentVolumeFraction>" << mQuiescentVolumeFraction << "</QuiescentVolumeFraction>\n";
    *rParamsFile << "\t\t\t<EquilibriumVolume>" << mEquilibriumVolume << "</EquilibriumVolume>\n";

    *rParamsFile << "\t\t\t<LabelledProbability>" << mLabelledProbability << "</LabelledProbability>\n";


    // Call method on direct parent class
    AbstractSimpleGenerationalCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptCellCycleModel)
