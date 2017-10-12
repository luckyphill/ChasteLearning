#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "FixedSequenceCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "VoronoiDataWriter.hpp"
#include "FakePetscSetup.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

class TestRunningMeshBasedSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    void TestMonolayer() throw(Exception)
    {
        CylindricalHoneycombMeshGenerator generator(5, 5);    // Parameters are: cells across, cells up
        //MutableMesh<2,2>* p_mesh = generator.GetMesh();
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<CellPtr> cells; 
        MAKE_PTR(TransitCellProliferativeType, p_transit_type); //Don't know why smart pointers need to be used here
        CellsGenerator<FixedSequenceCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());//, p_transit_type);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        CellCycleTimesGenerator::Instance()->GenerateCellCycleTimeSequence(); //needed for FixedSequenceCellCycleModel

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("GonnaStepOnTheGas");
        simulator.SetEndTime(12.0);

        simulator.SetSamplingTimestepMultiple(12);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();

        CellCycleTimesGenerator::Destroy(); //needed for FixedSequenceCellCycleModel

        //delete p_mesh;

        //TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);
        //TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 10.0, 1e-10);
    }

    void TestMonolayerWithGhostNodes() throw(Exception)
    {
        HoneycombMeshGenerator generator(4, 4, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices); //**Changed**//

        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TonightWellFly");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(12.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();

        //TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);
        //TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 10.0, 1e-10);
    }
};
