#include "DifferentiatedMembraneState.hpp"

DifferentiatedMembraneState::DifferentiatedMembraneState()
    : AbstractCellMutationState(0)
{}


#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(DifferentiatedMembraneState)
