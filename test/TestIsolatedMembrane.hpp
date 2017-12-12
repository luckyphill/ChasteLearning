#include <cxxtest/TestSuite.h> //Needed for all test files
#include "CellBasedSimulationArchiver.hpp" //Needed if we would like to save/load simulations
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing

/* The next set of classes are needed specifically for the simulation, which can be found in the core code. */

#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "FakePetscSetup.hpp"
#include "GeneralisedLinearSpringForce.hpp" //give a force to use between cells
#include "MembraneCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

#include "EpithelialLayerBasementMembraneForce.hpp"
#include "EpithelialLayerBasementMembraneForceModified.hpp"
#include "EpithelialLayerLinearSpringForce.hpp"

#include "LinearSpringForceMembraneCell.hpp"
#include "MembraneCellForce.hpp" // A force to restore the membrane to it's preferred shape
#include "NoCellCycleModel.hpp"

#include "CryptBoundaryCondition.hpp"
#include "BoundaryCellProperty.hpp"


class TestIsolatedMembrane : public AbstractCellBasedTestSuite
{
	public:
	void TestIsolatedFlatMembrane() throw(Exception)
	{
		unsigned cells_up = 10;
		unsigned cells_across = 30;

		double dt = 0.005;
		double end_time = 100;
		double sampling_multiple = 10;

		//Set all the spring stiffness variables
		double epithelialStiffness = 15.0; //Epithelial-epithelial spring connections
		double membraneStiffness = 60.0; //Stiffness of membrane to membrane spring connections
		double stromalStiffness = 15.0;

		double epithelialMembraneStiffness = 15.0; //Epithelial-non-epithelial spring connections
		double membraneStromalStiffness = 5.0; //Non-epithelial-non-epithelial spring connections
		double stromalEpithelialStiffness = 10.0;

		double torsional_stiffness = 5.0;
		double stiffness_ratio = 4.5; // For paneth cells
		
		double targetCurvatureStemStem = 0.2;
		double targetCurvatureStemTrans = 0; // Not implemented properly, so keep it the same as TransTrans for now
		double targetCurvatureTransTrans = 0;


		HoneycombMeshGenerator generator(cells_across, cells_up);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();

		std::vector<unsigned> initial_real_indices = generator.GetCellLocationIndices();
		std::vector<unsigned> real_indices;

		for (unsigned i = 0; i < initial_real_indices.size(); i++)
		{
			unsigned cell_index = initial_real_indices[i];
			double x = p_mesh->GetNode(cell_index)->rGetLocation()[0];
			double y = p_mesh->GetNode(cell_index)->rGetLocation()[1];

			//Make the curved crypt base
			if ( y > int(cells_up/2) && y < int(cells_up/2 + 1) && x <cells_across - 5 && x >5)
			{
				real_indices.push_back(cell_index);
			}
		}

		boost::shared_ptr<AbstractCellProperty> p_membrane = CellPropertyRegistry::Instance()->Get<MembraneCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

		std::vector<CellPtr> cells;

		for (unsigned i = 0; i<real_indices.size(); i++)
		{
			NoCellCycleModel* p_cycle_model = new NoCellCycleModel();
			CellPtr p_cell(new Cell(p_state, p_cycle_model));

			p_cell->SetCellProliferativeType(p_membrane);

			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell); 
		}

		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		OffLatticeSimulation<2> simulator(cell_population);

		simulator.SetOutputDirectory("TestIsolatedFlateMembrane");
        simulator.SetEndTime(end_time);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_multiple);

        MAKE_PTR(LinearSpringForceMembraneCell<2>, p_spring_force);
		p_spring_force->SetCutOffLength(1.5);
		//Set the spring stiffnesses
		p_spring_force->SetEpithelialSpringStiffness(epithelialStiffness);
		p_spring_force->SetMembraneSpringStiffness(membraneStiffness);
		p_spring_force->SetStromalSpringStiffness(stromalStiffness);
		p_spring_force->SetEpithelialMembraneSpringStiffness(epithelialMembraneStiffness);
		p_spring_force->SetMembraneStromalSpringStiffness(membraneStromalStiffness);
		p_spring_force->SetStromalEpithelialSpringStiffness(stromalEpithelialStiffness);

		p_spring_force->SetPanethCellStiffnessRatio(stiffness_ratio);
		simulator.AddForce(p_spring_force);

		MAKE_PTR(MembraneCellForce, p_membrane_force);
        p_membrane_force->SetBasementMembraneTorsionalStiffness(torsional_stiffness);
        p_membrane_force->SetTargetCurvatures(targetCurvatureStemStem, targetCurvatureStemTrans, targetCurvatureTransTrans);
        simulator.AddForce(p_membrane_force);

        simulator.Solve();

	}

};