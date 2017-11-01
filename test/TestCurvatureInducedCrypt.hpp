/* CREATED BY: Phillip Brown
 * A Chaste test that attempts to replicate the results by SJD
 */

#include <cxxtest/TestSuite.h> //Needed for all test files
#include "CellBasedSimulationArchiver.hpp" //Needed if we would like to save/load simulations
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing

/* The next set of classes are needed specifically for the simulation, which can be found in the core code. */

#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "DifferentiatedCellProliferativeType.hpp" //Stops cells from proliferating
#include "TransitCellProliferativeType.hpp"
#include "FakePetscSetup.hpp"
#include "GeneralisedLinearSpringForce.hpp" //give a force to use between cells
#include "UniformCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "EpithelialLayerAnoikisCellKiller.hpp"
#include "EpithelialLayerBasementMembraneForce.hpp"
#include "EpithelialLayerLinearSpringForce.hpp"
#include "TransitCellAnoikisResistantMutationState.hpp"
#include "DifferentiatedMembraneState.hpp" //not a very good name, supposed to be quick and dirty way to create a "membrane cell" by mutating a differentiated cell
#include "MembraneCellForce.hpp" // A force to restore the membrane to it's preferred shape


class TestCurvatureInducedCrypt : public AbstractCellBasedTestSuite
{
	public:
	void TestStartFromFlat() throw(Exception)
	{
		unsigned cells_up = 20;
		unsigned cells_across = 40;
		unsigned ghosts = 4;

		double dt = 0.005;
		double end_time = 10;
		double sampling_multiple = 10;

		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 10.0; //Epithelial-epithelial spring connections
		double epithelial_nonepithelial_stiffness = 5.0; //Epithelial-non-epithelial spring connections
		double nonepithelial_nonepithelial_stiffness = 10.0; //Non-epithelial-non-epithelial spring connections
		double membrane_stiffness = 10.0; //Stiffnes of mebrane to membrane spring connections
		double torsional_stiffness = 25.0;
		double stiffness_ratio = 4.5; // For paneth cells
		
		double targetCurvatureStemStem = 1/5;
		double targetCurvatureStemTrans = 1/10;
		double targetCurvatureTransTrans = 0;

		HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();

		//Sort through the indices and decide which ones are ghost nodes
		std::vector<unsigned> real_indices = generator.GetCellLocationIndices();

		boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
		boost::shared_ptr<AbstractCellProperty> p_trans_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();

		boost::shared_ptr<AbstractCellProperty> p_membrane_mutated = CellPropertyRegistry::Instance()->Get<DifferentiatedMembraneState>();

		//create a vector to store the cells; it is currently empty
		std::vector<CellPtr> cells;

		//go through the real indices and build some cells
		for (unsigned i = 0; i<real_indices.size(); i++)
		{	
			//Set cell cycle
			UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel();
			//p_cycle_model->SetCellCycleDuration(); //randomly chooses a duration
			double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //Randomly set birth time to stop pulsing behaviour
			p_cycle_model->SetBirthTime(-birth_time);
			

			CellPtr p_cell(new Cell(p_state, p_cycle_model));

			unsigned cell_index = real_indices[i];
			double x = p_mesh->GetNode(cell_index)->rGetLocation()[0];
			double y = p_mesh->GetNode(cell_index)->rGetLocation()[1];

			p_cell->SetCellProliferativeType(p_diff_type); //set the type to differentiated if it's not a ghost node - types will be reset as follows

			//add stems cells to the base of the crypt
			if ( y >= (cells_up - 1) * sqrt(3) /2 )
			{
				if (x > cells_across/2 - 5 && x < cells_across/2 + 5 && x>5 && x<cells_across - 5)
				{
					p_cell->SetCellProliferativeType(p_stem_type); //set the cell to stem if it's in the middle bunch
				} else 
				{
					p_cell->SetCellProliferativeType(p_trans_type);
				}
				
			}
			if( y < (cells_up - 1) * sqrt(3) /2  && y >= (cells_up - 2) * sqrt(3) /2)
			{
				p_cell->SetMutationState(p_membrane_mutated);
			}

			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		//Pull it all together
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		/* Define the simulation class. */
		OffLatticeSimulation<2> simulator(cell_population);

		MAKE_PTR(EpithelialLayerLinearSpringForce<2>, p_spring_force);
		p_spring_force->SetCutOffLength(1.5);
		//Set the spring stiffnesses
		p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness);
		p_spring_force->SetEpithelialNonepithelialSpringStiffness(epithelial_nonepithelial_stiffness);
		p_spring_force->SetNonepithelialNonepithelialSpringStiffness(nonepithelial_nonepithelial_stiffness);
		p_spring_force->SetMembraneSpringStiffness(membrane_stiffness);
		p_spring_force->SetPanethCellStiffnessRatio(stiffness_ratio);
		simulator.AddForce(p_spring_force);

		MAKE_PTR_ARGS(EpithelialLayerAnoikisCellKiller, p_anoikis_killer, (&cell_population));
		simulator.AddCellKiller(p_anoikis_killer);

		MAKE_PTR(MembraneCellForce, p_membrane_force);
        p_membrane_force->SetBasementMembraneTorsionalStiffness(torsional_stiffness);
        p_membrane_force->SetTargetCurvatures(targetCurvatureStemStem, targetCurvatureStemTrans, targetCurvatureTransTrans);
        simulator.AddForce(p_membrane_force);

		//Set output directory
		simulator.SetOutputDirectory("TestCurvatureInducedCrypt");
        simulator.SetEndTime(end_time);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_multiple);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();


	}
};