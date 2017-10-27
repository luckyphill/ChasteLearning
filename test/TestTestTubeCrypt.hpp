//Playing around with generating cells.
//Attempting to manually generate cells


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



class TestBasicTestTubeCrypt : public AbstractCellBasedTestSuite
{
	public:
	void TestGeneratedCells() throw(Exception)
	{
		unsigned cells_up = 30;
		unsigned cells_across = 30;
		unsigned ghosts = 4;
		double lumen_left_edge = 6.5;
		double lumen_right_edge = 13.5;
		double lumen_bottom = 4.5;


		//Copied wholesale from Axel
		//used to make the base of the test tube shaped crypt
		c_vector<double,2> circle_centre;
		circle_centre(0) = cells_across/2;
		circle_centre(1) = 10;

		double circle_radius = 5; //Size of hole

		double ring_width = 0.9;

		double dt = 0.01;
		double end_time = 20;
		double sampling_multiple = 100;
		//Basement membrane force parameters
		double bm_force = 6.0;
		double target_curvature = .15;
		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 5.0; //Epithelial-epithelial spring connections
		double epithelial_nonepithelial_stiffness = 5.0; //Epithelial-non-epithelial spring connections
		double nonepithelial_nonepithelial_stiffness = 15.0; //Non-epithelial-non-epithelial spring connections
		double membrane_stiffness = 25.0; //Stiffnes of mebrane to membrane spring connections
		//Set the stiffness ratio for Paneth cells to stem cells. This is the
		double stiffness_ratio = 4.5;
		//Start off with a mesh
		HoneycombMeshGenerator generator(cells_up, cells_across, ghosts);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();

	


		//Sort through the indices and decide which ones are ghost nodes
		std::vector<unsigned> initial_real_indices = generator.GetCellLocationIndices();
		std::vector<unsigned> real_indices;
		std::set<unsigned> real_indices_set;

		//for each index
			//check if index falls in the crypt lumen
				//if it doesn't, add to list of non-ghost indices

		for (unsigned i = 0; i < initial_real_indices.size(); i++)
		{
			unsigned cell_index = initial_real_indices[i];
			double x = p_mesh->GetNode(cell_index)->rGetLocation()[0];
			double y = p_mesh->GetNode(cell_index)->rGetLocation()[1];

			//Make the curved crypt base
			if ( (pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) > pow(circle_radius,2)) && y<= circle_centre[1])
			{
				real_indices.push_back(cell_index);
				real_indices_set.insert(cell_index);

			}

			if ( ((x <= circle_centre[0] - circle_radius -.5) || (x >= circle_centre[0] + circle_radius + .5)) && (y > circle_centre[1]))
			{
				real_indices.push_back(cell_index);
				real_indices_set.insert(cell_index);
			}
		}
		//When instantiating a cell, we need to provide a mutation state
		boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
		//Also need proliferative types:
		//Transit
		boost::shared_ptr<AbstractCellProperty> p_trans_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		//Differentiated
		boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
		//Stem
		boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();

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
			if ((pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) > pow(circle_radius,2)) && (pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) < pow(circle_radius + ring_width,2)) && y<= circle_centre[1])
			{
				p_cell->SetCellProliferativeType(p_stem_type); //set the cell to stem if it's at the base of the lumen
			}

			//add transit cells to the left edge of lumen
			if ( ((x <= circle_centre[0] - circle_radius -.5) && (x >= circle_centre[0] - circle_radius -1)) && (y > circle_centre[1]))
			{
				p_cell->SetCellProliferativeType(p_trans_type);
			}
			//add transit cells to the right edge of lumen
			if ( ((x >= circle_centre[0] + circle_radius +.5) && (x <= circle_centre[0] + circle_radius +1)) && (y > circle_centre[1]))
			{
				p_cell->SetCellProliferativeType(p_trans_type); 
			}
			

			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		//Pull it all together
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		/* Define the simulation class. */
		OffLatticeSimulation<2> simulator(cell_population);

		//Set output directory
		simulator.SetOutputDirectory("TestBasicTestTubeCrypt");
        simulator.SetEndTime(end_time);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_multiple);

        /* Add an anoikis-based cell killer. */
		MAKE_PTR_ARGS(EpithelialLayerAnoikisCellKiller, p_anoikis_killer, (&cell_population));
		simulator.AddCellKiller(p_anoikis_killer);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        //Basement membrane force
  //       MAKE_PTR(EpithelialLayerBasementMembraneForce, p_bm_force);
		// p_bm_force->SetBasementMembraneParameter(bm_force); //Equivalent to beta in SJD's papers
		// p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
		// simulator.AddForce(p_bm_force);


		MAKE_PTR(EpithelialLayerLinearSpringForce<2>, p_spring_force);
		p_spring_force->SetCutOffLength(1.5);
		//Set the spring stiffnesses
		p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness);
		p_spring_force->SetEpithelialNonepithelialSpringStiffness(epithelial_nonepithelial_stiffness);
		p_spring_force->SetNonepithelialNonepithelialSpringStiffness(nonepithelial_nonepithelial_stiffness);
		p_spring_force->SetMembraneSpringStiffness(membrane_stiffness);
		p_spring_force->SetPanethCellStiffnessRatio(stiffness_ratio);
		simulator.AddForce(p_spring_force);

		//mutate a cell so it does not die from anoikis
        boost::shared_ptr<AbstractCellProperty> p_state_mutated = CellPropertyRegistry::Instance()->Get<TransitCellAnoikisResistantMutationState>();
        //"mutate" a differentiated cell if it is under the monolayer
        boost::shared_ptr<AbstractCellProperty> p_membrane_mutated = CellPropertyRegistry::Instance()->Get<DifferentiatedMembraneState>();

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);

            if (node_index == 822) // Chosen from looking at the results from steady state
            {
                cell_iter->SetMutationState(p_state_mutated);
            }
            //if cell type is epithelial or stem, then make sure it's a monolayer and make all the connected differentiated cells a mutated type (if they haven't already)
            if (cell_iter->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() || cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
            {
	            std::set<unsigned> neighbouring_node_indices = cell_population.GetNeighbouringNodeIndices(node_index);
	            //loop to make sure we have a monolayer
	            unsigned real_neighbour_count=0;
	            for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
	         			iter != neighbouring_node_indices.end();
	         				++iter)
	    		{
	    			//count the number of non-ghost neighbours
	    			if (real_indices_set.find(*iter) != real_indices_set.end())
	    			{
	    				real_neighbour_count +=1;
	    			}
	    		}
	    		//if the cell doesn't have any ghost neighbours, then it's not a monolayer, so convert back to differntiated type
	    		if (real_neighbour_count == neighbouring_node_indices.size())
	    		{
	    			//std::cout <<"Changed "<< node_index << std::endl;
	    			cell_iter->SetCellProliferativeType(p_diff_type);
	    		} else {
	    		//loop to make the membrane cells
	    			//loop through neighbours
		    		for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
		         			iter != neighbouring_node_indices.end();
		         				++iter)
		    		{
		    			if (real_indices_set.find(*iter) != real_indices_set.end()) //make sure the node is not a ghost first
		    			{
		    				CellPtr neighbour = cell_population.GetCellUsingLocationIndex(*iter);
		    				//check if the cell type is differentiated, then if it is, add the "mutation"
			    			if (neighbour->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
			    			{
			    				//add mutation
			    				neighbour->SetMutationState(p_membrane_mutated);
			    			}
		    			}
		    		}
		    	}
	    	}
        }

        simulator.Solve();

	};
};