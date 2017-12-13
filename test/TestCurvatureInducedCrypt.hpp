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
#include "CylindricalHoneycombMeshGenerator.hpp" //Generates mesh
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
#include "EpithelialLayerBasementMembraneForce.hpp"
#include "EpithelialLayerBasementMembraneForceModified.hpp"
#include "EpithelialLayerLinearSpringForce.hpp"
#include "EpithelialLayerAnoikisCellKiller.hpp"

#include "AnoikisCellKillerMembraneCell.hpp"
#include "LinearSpringForceMembraneCell.hpp"
#include "MembraneCellForce.hpp" // A force to restore the membrane to it's preferred shape
#include "NoCellCycleModel.hpp"
#include "CryptBoundaryCondition.hpp"
#include "BoundaryCellProperty.hpp"


class TestCurvatureInducedCrypt : public AbstractCellBasedTestSuite
{
	public:
	void TestStartFromFlatMembraneCell() throw(Exception)
	{
		unsigned cells_up = 20;
		unsigned cells_across = 40;
		unsigned ghosts = 4;

		double dt = 0.001;
		double end_time = 10;
		double sampling_multiple = 1;

		//Set all the spring stiffness variables
		double epithelialStiffness = 15.0; //Epithelial-epithelial spring connections
		double membraneStiffness = 20.0; //Stiffness of membrane to membrane spring connections
		double stromalStiffness = 15.0;

		double epithelialMembraneStiffness = 10.0; //Epithelial-non-epithelial spring connections
		double membraneStromalStiffness = 15.0; //Non-epithelial-non-epithelial spring connections
		double stromalEpithelialStiffness = 10.0;

		double torsional_stiffness = 25.0;
		double stiffness_ratio = 4.5; // For paneth cells
		
		double targetCurvatureStemStem = 0.3;
		double targetCurvatureStemTrans = 0; // Not implemented properly, so keep it the same as TransTrans for now
		double targetCurvatureTransTrans = 0;

		CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
		// HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		// MutableMesh<2,2>* p_mesh = generator.GetMesh();

		//Sort through the indices and decide which ones are ghost nodes
		std::vector<unsigned> real_indices = generator.GetCellLocationIndices();

		boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
		boost::shared_ptr<AbstractCellProperty> p_trans_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_membrane = CellPropertyRegistry::Instance()->Get<MembraneCellProliferativeType>();

		

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

			if ( y >= (cells_up - 1.5) * sqrt(3) /2 )
			{
				if (x > cells_across/2 - 5 && x < cells_across/2 + 5)
				{
					p_cell->SetCellProliferativeType(p_stem_type); //set the cell to stem if it's in the middle bunch
				} else 
				{
					p_cell->SetCellProliferativeType(p_trans_type);
				}
				
			}
			if( y < (cells_up - 1.5) * sqrt(3) /2  && y >= (cells_up - 2) * sqrt(3) /2)
			{
				NoCellCycleModel* p_no_cycle_model = new NoCellCycleModel();
				p_cell->SetCellProliferativeType(p_membrane);
				p_cell->SetCellCycleModel(p_no_cycle_model);
			}

			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		//Pull it all together
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		/* Define the simulation class. */
		OffLatticeSimulation<2> simulator(cell_population);

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

		MAKE_PTR_ARGS(AnoikisCellKillerMembraneCell, p_anoikis_killer, (&cell_population));
		simulator.AddCellKiller(p_anoikis_killer);

		MAKE_PTR(MembraneCellForce, p_membrane_force);
        p_membrane_force->SetBasementMembraneTorsionalStiffness(torsional_stiffness);
        p_membrane_force->SetTargetCurvatures(targetCurvatureStemStem, targetCurvatureStemTrans, targetCurvatureTransTrans);
        simulator.AddForce(p_membrane_force);

		//Set output directory
		simulator.SetOutputDirectory("TestCurvatureInducedCryptMembraneCell");
        simulator.SetEndTime(end_time);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_multiple);

        simulator.Solve();

	}

	void xTestStartFromFlatMembraneForce() throw(Exception)
	{
		unsigned cells_up = 20;
		unsigned cells_across = 40;
		unsigned ghosts = 4;

		double dt = 0.001;
		double end_time = 100;
		double sampling_multiple = 100;

		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 15.0; //Epithelial-epithelial spring connections
		double epithelial_nonepithelial_stiffness = 10.0; //Epithelial-non-epithelial spring connections
		double nonepithelial_nonepithelial_stiffness = 15.0; //Non-epithelial-non-epithelial spring connections

		double stiffness_ratio = 4.5; // For paneth cells

		double bm_force = 10.0; //Set the basement membrane stiffness
		double target_curvature = 0.2; //Set the target curvature, i.e. how circular the layer wants to be

		std::cout << "1" << std::endl;
		CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
		// HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		// MutableMesh<2,2>* p_mesh = generator.GetMesh();

		//Sort through the indices and decide which ones are ghost nodes
		std::vector<unsigned> real_indices = generator.GetCellLocationIndices();
		std::cout << "2" << std::endl;
		boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
		boost::shared_ptr<AbstractCellProperty> p_trans_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();


		//create a vector to store the cells; it is currently empty
		std::vector<CellPtr> cells;
		std::cout << "3" << std::endl;
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
			if ( y >= (cells_up - 1.5) * sqrt(3) /2 )
			{
				if (x > cells_across/2 - 5 && x < cells_across/2 + 5)
				{
					p_cell->SetCellProliferativeType(p_stem_type); //set the cell to stem if it's in the middle bunch
				} else 
				{
					p_cell->SetCellProliferativeType(p_trans_type);
				}
				
			}

			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		std::cout << "4" << std::endl;
		//Pull it all together
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);
		std::cout << "5" << std::endl;
		cell_population.AddPopulationWriter<VoronoiDataWriter>();
		std::cout << "6" << std::endl;
		/* Define the simulation class. */
		OffLatticeSimulation<2> simulator(cell_population);
		std::cout << "7" << std::endl;
		MAKE_PTR(EpithelialLayerLinearSpringForce<2>, p_spring_force);
		p_spring_force->SetCutOffLength(1.5);
		//Set the spring stiffnesses
		p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness);
		p_spring_force->SetEpithelialNonepithelialSpringStiffness(epithelial_nonepithelial_stiffness);
		p_spring_force->SetNonepithelialNonepithelialSpringStiffness(nonepithelial_nonepithelial_stiffness);
		p_spring_force->SetPanethCellStiffnessRatio(stiffness_ratio);
		simulator.AddForce(p_spring_force);
		std::cout << "8" << std::endl;
		// Basement membrane force
        MAKE_PTR(EpithelialLayerBasementMembraneForceModified, p_bm_force);
		p_bm_force->SetBasementMembraneParameter(bm_force); //Equivalent to beta in SJD's papers
		p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
		simulator.AddForce(p_bm_force);
		
		std::cout << "9" << std::endl;
		MAKE_PTR_ARGS(EpithelialLayerAnoikisCellKiller, p_anoikis_killer, (&cell_population));
		simulator.AddCellKiller(p_anoikis_killer);
		std::cout << "10" << std::endl;

		// Stop it launching into the stratosphere
		MAKE_PTR_ARGS(CryptBoundaryCondition, p_bc, (&cell_population));
		simulator.AddCellPopulationBoundaryCondition(p_bc);

		//Set output directory
		simulator.SetOutputDirectory("TestCurvatureInducedCryptMembraneForce");
        simulator.SetEndTime(end_time);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_multiple);

        simulator.Solve();

	}

	void xTestIfNewMembraneForceWorks() throw(Exception)
	{
		// Build up a mesh then see how the forces get added in
		unsigned cells_up = 20;
		unsigned cells_across = 40;
		unsigned ghosts = 4;

		double dt = 0.001;
		double end_time = 100;
		double sampling_multiple = 100;

		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 15.0; //Epithelial-epithelial spring connections
		double epithelial_nonepithelial_stiffness = 15.0; //Epithelial-non-epithelial spring connections
		double nonepithelial_nonepithelial_stiffness = 15.0; //Non-epithelial-non-epithelial spring connections

		double stiffness_ratio = 4.5; // For paneth cells

		double bm_force = 10.0; //Set the basement membrane stiffness
		double target_curvature = 0.2; //Set the target curvature, i.e. how circular the layer wants to be

		std::cout << "1" << std::endl;
		CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
		// HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		// MutableMesh<2,2>* p_mesh = generator.GetMesh();

		//Sort through the indices and decide which ones are ghost nodes
		std::vector<unsigned> real_indices = generator.GetCellLocationIndices();
		std::cout << "2" << std::endl;
		boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
		boost::shared_ptr<AbstractCellProperty> p_trans_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();


		//create a vector to store the cells; it is currently empty
		std::vector<CellPtr> cells;
		std::cout << "3" << std::endl;
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
			if ( y >= (cells_up - 1.5) * sqrt(3) /2 )
			{
				if (x > cells_across/2 - 5 && x < cells_across/2 + 5)
				{
					p_cell->SetCellProliferativeType(p_stem_type); //set the cell to stem if it's in the middle bunch
				} else 
				{
					p_cell->SetCellProliferativeType(p_trans_type);
				}
				
			}

			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		std::cout << "4" << std::endl;
		//Pull it all together
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);
		std::cout << "5" << std::endl;
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		MAKE_PTR(EpithelialLayerBasementMembraneForceModified, p_bm_force);
		p_bm_force->SetBasementMembraneParameter(bm_force); //Equivalent to beta in SJD's papers
		p_bm_force->SetTargetCurvature(target_curvature);

		std::cout << "6" << std::endl;
		std::vector<c_vector<unsigned, 2> > pairs = p_bm_force->GetEpithelialStromalPairs(cell_population);
		
		for (unsigned i=0; i<pairs.size(); i++)
		{
			c_vector<double, 2> pair = pairs[i];
			std::cout << "Epithelial node: " << pair[0] << " Basement node: " << pair[1] << std::endl;
		}

		MAKE_PTR(EpithelialLayerBasementMembraneForce, p_bm_force2);
		p_bm_force2->SetBasementMembraneParameter(bm_force); //Equivalent to beta in SJD's papers
		p_bm_force2->SetTargetCurvature(target_curvature);
		std::vector<c_vector<unsigned, 2> > pairs2 = p_bm_force2->GetEpithelialGelPairs(cell_population);

		std::cout << "Old way \n" << std::endl;
		for (unsigned i=0; i<pairs2.size(); i++)
		{
			c_vector<double, 2> pair = pairs2[i];
			std::cout << "Epithelial node: " << pair[0] << " Basement node: " << pair[1] << std::endl;
		}
	}
};



