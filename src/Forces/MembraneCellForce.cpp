#include "MembraneCellForce.hpp"
#include "AbstractCellProperty.hpp"
#include "Debug.hpp"

/*
 * Created by: PHILLIP BROWN, 27/10/2017
 * Initial Structure borrows heavily from "EpithelialLayerBasementMembraneForce.cpp"
 * as found in the Chaste Paper Tutorials for the CryptFissionPlos2016 project
 */

/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
MembraneCellForce::MembraneCellForce()
   :  AbstractForce<2>(),
   mBasementMembraneTorsionalStiffness(DOUBLE_UNSET),
   mTargetCurvatureStemStem(DOUBLE_UNSET),
   mTargetCurvatureStemTrans(DOUBLE_UNSET),
   mTargetCurvatureTransTrans(DOUBLE_UNSET)
{
}

MembraneCellForce::~MembraneCellForce()
{

}

void MembraneCellForce::SetMembraneIndices(std::vector<unsigned> membraneIndices)
{
	mMembraneIndices = membraneIndices;
}

void MembraneCellForce::SetBasementMembraneTorsionalStiffness(double basementMembraneTorsionalStiffness)
{
	mBasementMembraneTorsionalStiffness = basementMembraneTorsionalStiffness;
}

double MembraneCellForce::GetBasementMembraneTorsionalStiffness()
{
	return mBasementMembraneTorsionalStiffness;
}


void MembraneCellForce::SetTargetCurvatures(double targetCurvatureStemStem, double targetCurvatureStemTrans, double targetCurvatureTransTrans)
{
	mTargetCurvatureStemStem = targetCurvatureStemStem;
	mTargetCurvatureStemTrans = targetCurvatureStemTrans;
	mTargetCurvatureTransTrans = targetCurvatureTransTrans;
}


double MembraneCellForce::GetTargetCurvatures(bool stem, bool trans)
{	
	//for a cell pair, enter true or false if either of them are attached to stem cells or transit cells
	//to determine appropriate target curvature
	assert((stem || trans)); //one of the options must be true
	if (stem && !trans)
	{
		return mTargetCurvatureStemStem;
	}
	if (!stem && trans)
	{
		return mTargetCurvatureTransTrans;
	}
	if (stem && trans)
	{
		return mTargetCurvatureStemTrans;
	}
	return 0;
}


/*
 * A method to find all the pairs of connections between healthy epithelial cells and labelled gel cells.
 * Returns a vector of node pairings, without repeats. The first of each pair is the epithelial node index,
 * and the second is the gel node index. Updating so that it also returns mutant-labelled cell pairs.
 */


double MembraneCellForce::GetAngleFromTriplet(AbstractCellPopulation<2>& rCellPopulation,
															c_vector<double, 2> leftNode,
															c_vector<double, 2> centreNode,
															c_vector<double, 2> rightNode)
{
	// Given three node which we know are neighbours, determine the angle their centres make
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	c_vector<double, 2> vector_AB = p_tissue->rGetMesh().GetVectorFromAtoB(leftNode,centreNode);
	c_vector<double, 2> vector_AC = p_tissue->rGetMesh().GetVectorFromAtoB(centreNode,rightNode);

	double inner_product_AB_AC = vector_AB[0] * vector_AC[0] + vector_AB[1] * vector_AC[1];

	double length_AB = norm_2(vector_AB);
	double length_AC = norm_2(vector_AC);

	double angle = acos(inner_product_AB_AC / (length_AB * length_AC));
	return angle;
	// Need to orient the vectors with respect to the lumen so we know which direction
	// THis is not done here, so must be done after
}

/*
* Function to return the curvature between three points parametrically - the midpoints of the springs connecting the
* transit cells to the differentiated cells. NB. The input arguments need to be in order from either left to right
* or right to left. If they are wrongly arranged (eg. middle, left, right) then you get a different curvature,
* but left->right = -(right-> left).
*/

double MembraneCellForce::FindParametricCurvature(AbstractCellPopulation<2>& rCellPopulation,
															c_vector<double, 2> leftMidpoint,
															c_vector<double, 2> centreMidpoint,
															c_vector<double, 2> rightMidpoint)
{
	//Get the relevant vectors (all possible differences)
	c_vector<double, 2> left_to_centre = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftMidpoint, centreMidpoint);
	c_vector<double, 2> centre_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centreMidpoint, rightMidpoint);
	c_vector<double, 2> left_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftMidpoint, rightMidpoint);

	// Firstly find the parametric intervals
	double left_s = sqrt(pow(left_to_centre[0],2) + pow(left_to_centre[1],2));
	double right_s = sqrt(pow(centre_to_right[0],2) + pow(centre_to_right[1],2));

	double sum_intervals = left_s + right_s;

	//Calculate finite difference of first derivatives
	double x_prime = (left_to_right[0])/sum_intervals;
	double y_prime = (left_to_right[1])/sum_intervals;

	//Calculate finite difference of second derivatives
	double x_double_prime = 2*(left_s*centre_to_right[0] - right_s*left_to_centre[0])/(left_s*right_s*sum_intervals);
	double y_double_prime = 2*(left_s*centre_to_right[1] - right_s*left_to_centre[1])/(left_s*right_s*sum_intervals);

	//Calculate curvature using formula
	double curvature = (x_prime*y_double_prime - y_prime*x_double_prime)/pow((pow(x_prime,2) + pow(y_prime,2)),3/2);

	return curvature;
}



//Method overriding the virtual method for AbstractForce. The crux of what really needs to be done.
void MembraneCellForce::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	// Need to determine the restoring force on the membrane putting it back to it's preferred shape

	// We loop over the epithelial-gel node pairs to find the force acting on that
	// epithelial node, and the direction in which it acts
	for (unsigned i=0; i<mMembraneIndices.size()-2; i++)
	{
		unsigned left_node = mMembraneIndices[i];
		unsigned centre_node = mMembraneIndices[i+1];
		unsigned right_node = mMembraneIndices[i+2];

		CellPtr left_cell = p_tissue->GetCellUsingLocationIndex(left_node);
		CellPtr centre_cell = p_tissue->GetCellUsingLocationIndex(centre_node);
		CellPtr right_cell = p_tissue->GetCellUsingLocationIndex(right_node);

		c_vector<double, 2> left_location = p_tissue->GetLocationOfCellCentre(left_cell);
		c_vector<double, 2> right_location = p_tissue->GetLocationOfCellCentre(right_cell);
		c_vector<double, 2> centre_location = p_tissue->GetLocationOfCellCentre(centre_cell);

		double current_angle = GetAngleFromTriplet(rCellPopulation, left_location, centre_location, right_location);
		double target_angle = 0.0;

		double torque = basementMembraneTorsionalStiffness * abs(current_angle - target_angle);

		/*
		 * Get the direction of force applied to each node
		 * Make a force vector
		 * Apply force to each node
		 */



   		//CellPtr p_cell_epithelial = p_tissue->GetCellUsingLocationIndex(epithelial_node_index);
   		//assert(p_cell_epithelial->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() == false);

		//CellPtr p_cell_gel = p_tissue->GetCellUsingLocationIndex(gel_node_index);
		//assert(p_cell_gel->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() == true);

		//c_vector<double, 2> epithelial_location = rCellPopulation.GetNode(epithelial_node_index)->rGetLocation();
		//c_vector<double, 2> gel_location = rCellPopulation.GetNode(gel_node_index)->rGetLocation();

		// The force due to the basal lamina acts along the spring connecting the epithelial and gel nodes, G->E direction
		//c_vector<double, 2> curvature_force_direction = p_tissue->rGetMesh().GetVectorFromAtoB(gel_location, epithelial_location);

		//double distance_between_nodes = norm_2(curvature_force_direction);
		//assert(distance_between_nodes > 0);
		//assert(!isnan(distance_between_nodes));

		//curvature_force_direction /= distance_between_nodes;

		//double curvature = GetCurvatureFromNodePair(rCellPopulation, epithelial_node_index, gel_node_index);

		//double basement_membrane_parameter = GetBasementMembraneParameter();

		//c_vector<double, 2> force_due_to_basement_membrane = basement_membrane_parameter*curvature*curvature_force_direction;

		// Add the force due to the basal lamina to the forces acting on that epithelial node
		//rCellPopulation.GetNode(epithelial_node_index)->AddAppliedForceContribution(force_due_to_basement_membrane);
	}

}

void MembraneCellForce::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<BasementMembraneTorsionalStiffness>"<<  mBasementMembraneTorsionalStiffness << "</BasementMembraneTorsionalStiffness> \n";
	//*rParamsFile <<  "\t\t\t<TargetCurvature>"<< mTargetCurvature << "</TargetCurvature> \n";

	// Call direct parent class
	AbstractForce<2>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MembraneCellForce)
