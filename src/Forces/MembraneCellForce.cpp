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


double MembraneCellForce::GetTargetCurvature(bool stem, bool trans)
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

	c_vector<double, 2> vector_AB = p_tissue->rGetMesh().GetVectorFromAtoB(centreNode,leftNode);
	c_vector<double, 2> vector_AC = p_tissue->rGetMesh().GetVectorFromAtoB(centreNode,rightNode);

	double inner_product_AB_AC = vector_AB[0] * vector_AC[0] + vector_AB[1] * vector_AC[1];
	//std::cout << "Inner Product: " << inner_product_AB_AC << std::endl;
	double length_AB = norm_2(vector_AB);
	double length_AC = norm_2(vector_AC);
	//std::cout << "Length C-L: " << length_AB << "\nLength C-R: " << length_AC << std::endl;

	double acos_arg = inner_product_AB_AC / (length_AB * length_AC);

	//std::cout << "Argument: " << acos_arg << setprecision(17) << std::endl;

	double angle = acos(acos_arg);
	// Occasionally the argument steps out of the bounds for acos, for instance -1.0000000000000002
	// This is enough to make the acos function return nans
	// This line of code is not ideal, but catches the error for the time being
	if (isnan(angle) && acos_arg > -1.00000000000005) 
	{
		return acos(-1);
	}
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
															c_vector<double, 2> leftCell,
															c_vector<double, 2> centreCell,
															c_vector<double, 2> rightCell)
{
	//Get the relevant vectors (all possible differences)
	c_vector<double, 2> left_to_centre = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftCell, centreCell);
	c_vector<double, 2> centre_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centreCell, rightCell);
	c_vector<double, 2> left_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftCell, rightCell);

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


double MembraneCellForce::GetTargetAngle(AbstractCellPopulation<2>& rCellPopulation, CellPtr centre_cell,
																		c_vector<double, 2> leftCell,
																		c_vector<double, 2> centreCell,
																		c_vector<double, 2> rightCell)
{
	// Returns the angle that we're aiming for
	// At the moment, it doesn't handle membrane cells with both types separately, but treats them like  they're attached to transit cells
	MeshBasedCellPopulation<2>* cell_population = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	bool contact_with_stem = false;
	bool contact_with_trans = false;

	//std::cout << "About to assign centre cell index" <<std::endl;
	unsigned centre_cell_index = cell_population->GetLocationIndexUsingCell(centre_cell);
	//std::cout << "Successfully assigned centre cell index" <<std::endl;
	//std::cout << "About to get neighbours" <<std::endl;
	std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(centre_cell_index);
	//std::cout << "Successfully got neighbours" <<std::endl;

	//std::cout << "About to loop neighbours" <<std::endl;
	for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
	         			iter != neighbouring_node_indices.end();
	         				++iter)
	{
		if (!cell_population->IsGhostNode(*iter))
		{		
			CellPtr neighbour = cell_population->GetCellUsingLocationIndex(*iter);
			//check if the cell type is differentiated, then if it is, add the "mutation"
			if (neighbour->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
			{
				contact_with_trans = true;
			}
			if (neighbour->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
			{
				contact_with_stem = true;
			}
		}
	}
	//std::cout << "Successfully looped neighbours" <<std::endl;

	assert((contact_with_stem || contact_with_trans));

	//double target_curvature = GetTargetCurvature(contact_with_stem, contact_with_trans);

	double target_angle = M_PI; // Assume we're dealing with only transit cells for the now

	c_vector<double, 2> vector_AB = cell_population->rGetMesh().GetVectorFromAtoB(centreCell,leftCell);
	c_vector<double, 2> vector_AC = cell_population->rGetMesh().GetVectorFromAtoB(centreCell,rightCell);

	double length_AB = norm_2(vector_AB);
	double length_AC = norm_2(vector_AC);


	if (contact_with_stem && !contact_with_trans)
	{
		//target_angle = acos(length_AC * mTargetCurvatureStemStem / 2) + acos(length_AB * mTargetCurvatureStemStem / 2);
		target_angle = 3.0;
	}

	// if (contact_with_stem && contact_with_trans)
	// {
	// 	target_angle = acos(length_AC * target_curvature) + acos(length_AB * target_curvature);
	// }
	return target_angle;
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
		//std::cout << "\nNode: " << centre_node  <<std::endl;

		//std::cout << "About to assign cells" <<std::endl;
		CellPtr left_cell = p_tissue->GetCellUsingLocationIndex(left_node);
		CellPtr centre_cell = p_tissue->GetCellUsingLocationIndex(centre_node);
		CellPtr right_cell = p_tissue->GetCellUsingLocationIndex(right_node);
		//std::cout << "Successfully assigned cells" <<std::endl;

		c_vector<double, 2> left_location = p_tissue->GetLocationOfCellCentre(left_cell);
		c_vector<double, 2> right_location = p_tissue->GetLocationOfCellCentre(right_cell);
		c_vector<double, 2> centre_location = p_tissue->GetLocationOfCellCentre(centre_cell);

		
		//std::cout << "About to get angle and curvature" <<std::endl;
		double current_angle = GetAngleFromTriplet(rCellPopulation, left_location, centre_location, right_location);
		double current_curvature = FindParametricCurvature(rCellPopulation, left_location, centre_location, right_location);
		//std::cout << "Successfully got angle and curvature" << std::endl;
		
		if (abs(current_curvature) < 1e-5)
		{
			// Close enough
			current_curvature = 0.0;
			// We need to use the sign of the curvature to determine the angle correctly
			// Extrememly small curvatures due to precision errors might play havock with this
		}

		// The method of calculating the angle is not oriented by the lumen, so need to adjust
		if (current_curvature < 0)
		{
			if (centre_node == 507)
			{
				std::cout << "Changing current angle\nStarting: " << current_angle << std::endl;
			}
			current_angle = 2 * M_PI - current_angle;
			if (centre_node == 507)
			{
				std::cout << "Ending: " << current_angle << std::endl;
			}
		}

		//std::cout << "About to get target angle"  <<std::endl;
		double target_angle = GetTargetAngle(rCellPopulation, centre_cell, left_location, centre_location, right_location);
		//std::cout << "Successfully got target angle"  <<std::endl;

		//std::cout<< "Angle: " << current_angle << "\nTarget Angle: " << target_angle << std::endl;
		//std::cout << "Left Node: " << left_location[0] << ", " << left_location[1] << "\nCentre Node: " << centre_location[0]<< ", " << centre_location[1] << "\nRight Node: " << right_location[0]<< ", " << right_location[1] << std::endl;

		double torque = mBasementMembraneTorsionalStiffness * (current_angle - target_angle); // Positive torque means force points into lumen
		//std::cout << "Torque: " << torque << std::endl;
		c_vector<double, 2> vector_CL = p_tissue->rGetMesh().GetVectorFromAtoB(centre_location,left_location);
		c_vector<double, 2> vector_CR = p_tissue->rGetMesh().GetVectorFromAtoB(centre_location,right_location);
		c_vector<double, 2> vector_LR = p_tissue->rGetMesh().GetVectorFromAtoB(left_location,right_location); // Used for determining where lumen is



		double length_CL = norm_2(vector_CL);
		double length_CR = norm_2(vector_CR);

		// Determine the force vectors applied to the left and right nodes
		double forceMagnitudeLeft = length_CL * torque;
		double forceMagnitudeRight = length_CR * torque;
		//std::cout << "Force Magnitude: " << forceMagnitudeLeft << std::endl;

		c_vector<double, 2> forceDirectionLeft;
		c_vector<double, 2> forceDirectionRight;

		// Use the CL and CR vectors to determine the line that the force will act on
		// If we have a vector (a, b), then the vector (b, -a) is perpendicular and creates a clockwise rotation when added to the end of (a,b)
		// while (-b, a) creates an anticlockwise rotation
		// forceDirectionLeft will always end up pointing into the lumen, and forceDirectionRight will always point out
		// Given we have decided that the actual direction of the force is encoded in the sign on the torque this is all we need to do
		forceDirectionLeft[0] = vector_CL[1] / length_CL;
		forceDirectionLeft[1] = - vector_CL[0] / length_CL;

		forceDirectionRight[0] = - vector_CR[1] / length_CR;
		forceDirectionRight[1] = vector_CR[0] / length_CR;

		c_vector<double, 2> forceVectorLeft = forceMagnitudeLeft * forceDirectionLeft;
		c_vector<double, 2> forceVectorRight = forceMagnitudeRight * forceDirectionRight;

		if (centre_node == 507)
		{
			std::cout << "\nNode: " << centre_node  <<std::endl;
			std::cout << "Torque: " << torque << std::endl;
			std::cout<< "Angle: " << current_angle << "\nTarget Angle: " << target_angle << std::endl;
			std::cout<< "Curvature: " << current_curvature << std::endl;
			std::cout << "forceVectorLeft: " << forceVectorLeft[0]<< ", " << forceVectorLeft[1] <<std::endl;
			std::cout << "forceMagnitudeLeft: " << forceMagnitudeLeft <<std::endl;
			std::cout << "forceVectorRight: " << forceVectorRight[0]<< ", " << forceVectorRight[1] <<std::endl;
			std::cout << "forceMagnitudeRight: " << forceMagnitudeRight <<std::endl;

		}

		//std::cout << "About to add force contribution"  <<std::endl;
		rCellPopulation.GetNode(left_node)->AddAppliedForceContribution(forceVectorLeft);
		rCellPopulation.GetNode(right_node)->AddAppliedForceContribution(forceVectorRight);
		//std::cout << "Successfully added force contribution"  <<std::endl;
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
