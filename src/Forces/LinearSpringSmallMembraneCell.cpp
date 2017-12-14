/*
MODIFIED BY PHILLIP BROWN: 27/10/2017
- Added mutation that turns a differentiated cell into a "membrane cell" in 
in order to test a method of introducing a membrane
- The modifications here only change the way the "mutant" cells interact with each
other. Otherwise they are still considered "differentiated" cells for other interactions
MODIFICATIONS around lines 48, 211, 259, 410, 448
MODIFIED BY AXEL ALMET: 23/12/14
Copyright (c) 2005-2014, University of Oxford.
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

#include "LinearSpringSmallMembraneCell.hpp"
#include "IsNan.hpp"
#include "AbstractCellProperty.hpp"

#include "PanethCellMutationState.hpp"
#include "MembraneCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::LinearSpringSmallMembraneCell()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
    mEpithelialSpringStiffness(15.0), // Epithelial covers stem and transit
    mMembraneSpringStiffness(15.0),
    mStromalSpringStiffness(15.0), // Stromal is the differentiated "filler" cells
    mEpithelialMembraneSpringStiffness(15.0),
    mMembraneStromalSpringStiffness(15.0),
    mStromalEpithelialSpringStiffness(15.0),
    mEpithelialRestLength(1.0),
    mMembraneRestLength(1.0),
    mStromalRestLength(1.0),
    mEpithelialMembraneRestLength(1.0),
    mMembraneStromalRestLength(1.0),
    mStromalEpithelialRestLength(1.0),
    mEpithelialCutOffLength(1.5), // Epithelial covers stem and transit
    mMembraneCutOffLength(1.5),
    mStromalCutOffLength(1.5), // Stromal is the differentiated "filler" cells
    mEpithelialMembraneCutOffLength(1.5),
    mMembraneStromalCutOffLength(1.5),
    mStromalEpithelialCutOffLength(1.5),
    mPanethCellStiffnessRatio(1.0)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                                                     unsigned nodeBGlobalIndex,
                                                                                     AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                                                     bool isCloserThanRestLength)
{
    return 1.0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::~LinearSpringSmallMembraneCell()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<SPACE_DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<SPACE_DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, SPACE_DIM> node_a_location = p_node_a->rGetLocation();
    c_vector<double, SPACE_DIM> node_b_location = p_node_b->rGetLocation();


    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, SPACE_DIM> unitForceDirection;

    unitForceDirection = rCellPopulation.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unitForceDirection);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unitForceDirection /= distance_between_nodes;

    /*
     * Calculate the rest length of the spring connecting the two nodes with a default
     * value of 1.0.
     */

    // We have three types of cells, with 6 different possible pairings as demarked by the 6 different spring stiffnesses
    // Need to check which types we have and set spring_constant accordingly
    // There is also a method that gives the possibilty of a variable spring constant based on whether the spring is in tension or compression
    // this is not implemented here

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    // First, determine what we've got
    bool membraneA = p_cell_A->GetCellProliferativeType()->IsType<MembraneCellProliferativeType>();
    bool membraneB = p_cell_B->GetCellProliferativeType()->IsType<MembraneCellProliferativeType>();

    bool stromalA = p_cell_A->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>();
    bool stromalB = p_cell_B->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>();

    bool epiA = ( p_cell_A->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() || p_cell_A->GetCellProliferativeType()->IsType<StemCellProliferativeType>() );
    bool epiB = ( p_cell_B->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() || p_cell_B->GetCellProliferativeType()->IsType<StemCellProliferativeType>() );


    double rest_length_final = 1.0;
    rest_length_final = static_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)->GetRestLength(nodeAGlobalIndex, nodeBGlobalIndex);

    double spring_constant = 0.0;

    // Determine rest lengths and spring stiffnesses
    if (membraneA)
    {
        if (membraneB)
        {
            rest_length_final = mMembraneRestLength;
            spring_constant = mMembraneSpringStiffness;
        }
        if (stromalB)
        {   
            if (distance_between_nodes >= mMembraneStromalCutOffLength)
            {
                return zero_vector<double>(SPACE_DIM);
            }
            rest_length_final = mMembraneStromalRestLength;
            spring_constant = mMembraneStromalSpringStiffness;
        }
        if (epiB)
        {
            rest_length_final = mEpithelialMembraneRestLength;
            spring_constant = mEpithelialMembraneSpringStiffness;
        }
    }

    if (stromalA)
    {
        if (membraneB)
        {
            if (distance_between_nodes >= mMembraneStromalCutOffLength)
            {
                return zero_vector<double>(SPACE_DIM);
            }
            rest_length_final = mMembraneStromalRestLength;
            spring_constant = mMembraneStromalSpringStiffness;
        }
        if (stromalB)
        {
            rest_length_final = mStromalRestLength;
            spring_constant = mStromalSpringStiffness;
        }
        if (epiB)
        {
            rest_length_final = mStromalEpithelialRestLength;
            spring_constant = mStromalEpithelialSpringStiffness;
        }
    }

    if (epiA)
    {
        if (membraneB)
        {
            rest_length_final = mEpithelialMembraneRestLength;
            spring_constant = mEpithelialMembraneSpringStiffness;
        }
        if (stromalB)
        {
            rest_length_final = mStromalEpithelialRestLength;
            spring_constant = mStromalEpithelialSpringStiffness;
        }
        if (epiB)
        {
            rest_length_final = mEpithelialRestLength;
            spring_constant = mEpithelialSpringStiffness;
        }
    }

    assert(spring_constant > 0);
    double rest_length = rest_length_final;

    double ageA = p_cell_A->GetAge();
    double ageB = p_cell_B->GetAge();

    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));

    /*
     * If the cells are both newly divided, then the rest length of the spring
     * connecting them grows linearly with time, until 1 hour after division.
     */
    if (ageA < mMeinekeSpringGrowthDuration && ageB < mMeinekeSpringGrowthDuration)
    {
        AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

        std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

        if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
        {
            // Spring rest length increases from a small value to the normal rest length over 1 hour
            double lambda = mMeinekeDivisionRestingSpringLength;
            rest_length = lambda + (rest_length_final - lambda) * ageA/mMeinekeSpringGrowthDuration;
        }
        if (ageA + SimulationTime::Instance()->GetTimeStep() >= mMeinekeSpringGrowthDuration)
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

    double length_change = distance_between_nodes - rest_length;

    // std::cout << "Node pair: " << nodeAGlobalIndex << ", " << nodeBGlobalIndex << std::endl;
    // std::cout << "Spring constant: " << spring_constant << std::endl;
    // std::cout << "Direction x: " << unitForceDirection[0] << " Direction y: " << unitForceDirection[1] << std::endl;
    // std::cout << "Length change: " << length_change << std::endl;
    return spring_constant * length_change * unitForceDirection;

    

}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::GetEpithelialSpringStiffness()
{
    return mEpithelialSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::GetMembraneSpringStiffness()
{
    return mMembraneSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::GetStromalSpringStiffness()
{
    return mStromalSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::GetEpithelialMembraneSpringStiffness()
{
    return mEpithelialMembraneSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::GetMembraneStromalSpringStiffness()
{
    return mMembraneStromalSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::GetStromalEpithelialSpringStiffness()
{
    return mStromalEpithelialSpringStiffness;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::GetMeinekeDivisionRestingSpringLength()
{
    return mMeinekeDivisionRestingSpringLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::GetMeinekeSpringGrowthDuration()
{
    return mMeinekeSpringGrowthDuration;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::GetPanethCellStiffnessRatio()
{
    return mPanethCellStiffnessRatio;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetEpithelialSpringStiffness(double epithelialSpringStiffness)
{
    assert(epithelialSpringStiffness> 0.0);
    mEpithelialSpringStiffness = epithelialSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetMembraneSpringStiffness(double membraneSpringStiffness)
{
    assert(membraneSpringStiffness > 0.0);
    mMembraneSpringStiffness = membraneSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetStromalSpringStiffness(double stromalSpringStiffness)
{
    assert(stromalSpringStiffness > 0.0);
    mStromalSpringStiffness = stromalSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetEpithelialMembraneSpringStiffness(double epithelialMembraneSpringStiffness)
{
    assert(epithelialMembraneSpringStiffness > 0.0);
    mEpithelialMembraneSpringStiffness = epithelialMembraneSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetMembraneStromalSpringStiffness(double membraneStromalSpringStiffness)
{
    assert(membraneStromalSpringStiffness > 0.0);
    mMembraneStromalSpringStiffness = membraneStromalSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetStromalEpithelialSpringStiffness(double stromalEpithelialSpringStiffness)
{
    assert(stromalEpithelialSpringStiffness > 0.0);
    mStromalEpithelialSpringStiffness = stromalEpithelialSpringStiffness;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetEpithelialRestLength(double epithelialRestLength)
{
    assert(epithelialRestLength> 0.0);
    mEpithelialRestLength = epithelialRestLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetMembraneRestLength(double membraneRestLength)
{
    assert(membraneRestLength > 0.0);
    mMembraneRestLength = membraneRestLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetStromalRestLength(double stromalRestLength)
{
    assert(stromalRestLength > 0.0);
    mStromalRestLength = stromalRestLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetEpithelialMembraneRestLength(double epithelialMembraneRestLength)
{
    assert(epithelialMembraneRestLength > 0.0);
    mEpithelialMembraneRestLength = epithelialMembraneRestLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetMembraneStromalRestLength(double membraneStromalRestLength)
{
    assert(membraneStromalRestLength > 0.0);
    mMembraneStromalRestLength = membraneStromalRestLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetStromalEpithelialRestLength(double stromalEpithelialRestLength)
{
    assert(stromalEpithelialRestLength > 0.0);
    mStromalEpithelialRestLength = stromalEpithelialRestLength;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetEpithelialCutOffLength(double epithelialCutOffLength)
{
    assert(epithelialCutOffLength> 0.0);
    mEpithelialCutOffLength = epithelialCutOffLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetMembraneCutOffLength(double membraneCutOffLength)
{
    assert(membraneCutOffLength > 0.0);
    mMembraneCutOffLength = membraneCutOffLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetStromalCutOffLength(double stromalCutOffLength)
{
    assert(stromalCutOffLength > 0.0);
    mStromalCutOffLength = stromalCutOffLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetEpithelialMembraneCutOffLength(double epithelialMembraneCutOffLength)
{
    assert(epithelialMembraneCutOffLength > 0.0);
    mEpithelialMembraneCutOffLength = epithelialMembraneCutOffLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetMembraneStromalCutOffLength(double membraneStromalCutOffLength)
{
    assert(membraneStromalCutOffLength > 0.0);
    mMembraneStromalCutOffLength = membraneStromalCutOffLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetStromalEpithelialCutOffLength(double stromalEpithelialCutOffLength)
{
    assert(stromalEpithelialCutOffLength > 0.0);
    mStromalEpithelialCutOffLength = stromalEpithelialCutOffLength;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetMeinekeDivisionRestingSpringLength(double divisionRestingSpringLength)
{
    assert(divisionRestingSpringLength <= 1.0);
    assert(divisionRestingSpringLength >= 0.0);

    mMeinekeDivisionRestingSpringLength = divisionRestingSpringLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetMeinekeSpringGrowthDuration(double springGrowthDuration)
{
    assert(springGrowthDuration >= 0.0);

    mMeinekeSpringGrowthDuration = springGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::SetPanethCellStiffnessRatio(double panethCellStiffnessRatio)
{
    assert(panethCellStiffnessRatio >= 0.0);

    mPanethCellStiffnessRatio = panethCellStiffnessRatio;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringSmallMembraneCell<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<EpithelialSpringStiffness>" << mEpithelialSpringStiffness << "</EpithelialSpringStiffness>\n";
    *rParamsFile << "\t\t\t<MembraneSpringStiffness>" << mMembraneSpringStiffness << "</MembraneSpringStiffness>\n";
    *rParamsFile << "\t\t\t<StromalSpringStiffness>" << mStromalSpringStiffness << "</StromalSpringStiffness>\n";
    *rParamsFile << "\t\t\t<EpithelialMembraneSpringStiffness>" << mEpithelialMembraneSpringStiffness << "</EpithelialMembraneSpringStiffness>\n";
    *rParamsFile << "\t\t\t<MembranetromalSpringStiffness>" << mMembraneStromalSpringStiffness << "</MembranetromalSpringStiffness>\n";
    *rParamsFile << "\t\t\t<StromalEpithelialSpringStiffness>" << mStromalEpithelialSpringStiffness << "</StromalEpithelialSpringStiffness>\n";

    *rParamsFile << "\t\t\t<EpithelialRestLength>" << mEpithelialRestLength << "</EpithelialRestLength>\n";
    *rParamsFile << "\t\t\t<MembraneRestLength>" << mMembraneRestLength << "</MembraneRestLength>\n";
    *rParamsFile << "\t\t\t<StromalRestLength>" << mStromalRestLength << "</StromalRestLength>\n";
    *rParamsFile << "\t\t\t<EpithelialMembraneRestLength>" << mEpithelialMembraneRestLength << "</EpithelialMembraneRestLength>\n";
    *rParamsFile << "\t\t\t<MembranetromalRestLength>" << mMembraneStromalRestLength << "</MembranetromalRestLength>\n";
    *rParamsFile << "\t\t\t<StromalEpithelialRestLength>" << mStromalEpithelialRestLength << "</StromalEpithelialRestLength>\n";

    *rParamsFile << "\t\t\t<EpithelialCutOffLength>" << mEpithelialCutOffLength << "</EpithelialCutOffLength>\n";
    *rParamsFile << "\t\t\t<MembraneCutOffLength>" << mMembraneCutOffLength << "</MembraneCutOffLength>\n";
    *rParamsFile << "\t\t\t<StromalCutOffLength>" << mStromalCutOffLength << "</StromalCutOffLength>\n";
    *rParamsFile << "\t\t\t<EpithelialMembraneCutOffLength>" << mEpithelialMembraneCutOffLength << "</EpithelialMembraneCutOffLength>\n";
    *rParamsFile << "\t\t\t<MembranetromalCutOffLength>" << mMembraneStromalCutOffLength << "</MembranetromalCutOffLength>\n";
    *rParamsFile << "\t\t\t<StromalEpithelialCutOffLength>" << mStromalEpithelialCutOffLength << "</StromalEpithelialCutOffLength>\n";

    *rParamsFile << "\t\t\t<MeinekeDivisionRestingSpringLength>" << mMeinekeDivisionRestingSpringLength << "</MeinekeDivisionRestingSpringLength>\n";
    *rParamsFile << "\t\t\t<MeinekeSpringGrowthDuration>" << mMeinekeSpringGrowthDuration << "</MeinekeSpringGrowthDuration>\n";
    *rParamsFile << "\t\t\t<PanethCellStiffnessRatio>" << mPanethCellStiffnessRatio << "</PanethCellStiffnessRatio>\n";

    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class LinearSpringSmallMembraneCell<1,1>;
template class LinearSpringSmallMembraneCell<1,2>;
template class LinearSpringSmallMembraneCell<2,2>;
template class LinearSpringSmallMembraneCell<1,3>;
template class LinearSpringSmallMembraneCell<2,3>;
template class LinearSpringSmallMembraneCell<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(LinearSpringSmallMembraneCell)