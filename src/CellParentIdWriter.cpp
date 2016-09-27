/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "CellParentIdWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellParentId.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellParentIdWriter<ELEMENT_DIM, SPACE_DIM>::CellParentIdWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("parentcell.dat")
{
    this->mVtkCellDataName = "Cell Parent IDs";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellParentIdWriter<ELEMENT_DIM, SPACE_DIM>::GetCellParentId(CellPtr pCell)
{
    // CellPropertyCollection cell_parent_id_collection = pCell->rGetCellPropertyCollection().GetPropertiesType<CellParentId>();
    // assert(cell_parent_id_collection.GetSize() == 1);
    // boost::shared_ptr<CellParentId> p_cell_parent_id = boost::static_pointer_cast<CellParentId>(cell_parent_id_collection.GetProperty());
    // return p_cell_parent_id->GetParentId();
    return pCell->GetCellData()->GetItem("parent_id");
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> CellParentIdWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDivisionLocation(CellPtr pCell)
{
    c_vector<double, 3> location;
    
    location[0] = pCell->GetCellData()->GetItem("division_location_x");
    location[1] = pCell->GetCellData()->GetItem("division_location_y");
    location[2] = pCell->GetCellData()->GetItem("division_location_z");

    return location;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellParentIdWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double parent_cell_id = GetCellParentId(pCell);
    return parent_cell_id;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellParentIdWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned cell_id = pCell->GetCellId();
    int parent_cell_id = GetCellParentId(pCell);
    
    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    *this->mpOutStream << " " << cell_id << " " << parent_cell_id << " " << location_index;

    // Output coordinates of division
    c_vector<double, 3> coords = GetCellDivisionLocation(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << " " << coords[i];
    }

    *this->mpOutStream << " " << pCell->GetCellData()->GetItem("division_time");


}

// Explicit instantiation
template class CellParentIdWriter<1,1>;
template class CellParentIdWriter<1,2>;
template class CellParentIdWriter<2,2>;
template class CellParentIdWriter<1,3>;
template class CellParentIdWriter<2,3>;
template class CellParentIdWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellParentIdWriter)
