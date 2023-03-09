/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = custom_update_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = ecm_function; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 

	std::string mat_filename = parameters.strings("substrate_file");
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	read_custom_microenvironment_from_matlab(mat_filename);
	add_low_cut();

	return; 
}

void setup_tissue( void )
{
	// for(int n=0 ; n < microenvironment.number_of_voxels() ; n++){
	// 	std::cout << microenvironment(n) << std::endl;
	// }
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void ecm_function(Cell *pCell, Phenotype &phenotype, double dt) 
{
	// finds substrate indexes
	int ecm_1 = microenvironment.find_density_index("ECM_collagen");
	int ecm_2 = microenvironment.find_density_index("ECM_fibrin");
	int i = pCell -> get_current_mechanics_voxel_index();
	// finds concentrations in nearest voxel
	double conc_1 = pCell->get_microenvironment()->nearest_density_vector(i)[ecm_1];
	double conc_2 = pCell->get_microenvironment()->nearest_density_vector(i)[ecm_2];
	
	double r = UniformRandom();
	
	//checks for non diffusing ecm and gives cell a probability of degrading it
	if ((conc_1 >= parameters.doubles("ecm1_threshold")) && (r > parameters.doubles("secretion_probability")))
	{
		// cell uptakes non diffusing ecm and secretes diffusing ecm
		pCell->phenotype.secretion.uptake_rates[ecm_1] = parameters.doubles("ecm1_uptake");
		pCell->phenotype.secretion.secretion_rates[ecm_2] = parameters.doubles("ecm2_secretion");
	}
		
	else 
	{
		// otherwise no secretion or uptake
		pCell->phenotype.secretion.uptake_rates[ecm_1] = 0.0;
		pCell->phenotype.secretion.secretion_rates[ecm_2] = 0.0;
	}
}

double EPSILON = 1e-16;

/* Calculate repulsion/adhesion between agent and ecm according to its local density */
void add_ecm_interaction(Cell* pC, int index_ecm, int index_voxel )
{
	// Check if there is ECM material in given voxel
	//double dens2 = get_microenvironment()->density_vector(index_voxel)[index_ecm];
	double dens = pC->get_microenvironment()->nearest_density_vector(index_voxel)[index_ecm];
	double ecmrad = sqrt(3.0) * pC->get_microenvironment()->mesh.dx / 2.0;
	// if voxel is "full", density is 1
	dens = std::min( dens, 1.0 ); 
	if ( dens > EPSILON )
	{
		// Distance between agent center and ECM voxel center
		pC->displacement = pC->position - pC->get_container()->underlying_mesh.voxels[index_voxel].center;
		double distance = norm(pC->displacement);
		// Make sure that the distance is not zero
		distance = std::max(distance, EPSILON);
		
		double dd = pC->phenotype.geometry.radius + ecmrad;  
		double dnuc = pC->phenotype.geometry.nuclear_radius + ecmrad;  

		double tmp_r = 0;
		// Cell overlap with ECM node, add a repulsion term
		if ( distance < dd )
		{
			// repulsion stronger if nucleii overlap, see Macklin et al. 2012, 2.3.1
			if ( distance < dnuc )
			{
				double M = 1.0;
				double c = 1.0 - dnuc/dd;
				c *= c;
				c -= M;
				tmp_r = c*distance/dnuc + M;
			}
			else
			{
				tmp_r = ( 1 - distance / dd );
				tmp_r *= tmp_r;
			}
			tmp_r *= dens * PhysiCell::parameters.doubles("cell_ecm1_repulsion");
		}

		// Cell adherence to ECM through integrins
		double max_interactive_distance = (PhysiCell::parameters.doubles("max_interaction_factor")*pC->phenotype.geometry.radius) + ecmrad;
		if ( distance < max_interactive_distance ) 
		{	
			double temp_a = 1 - distance/max_interactive_distance; 
			temp_a *= temp_a; 
			/* \todo change dens with a maximal density ratio ? */
			// temp_a *= dens * ( static_cast<Cell*>(this) )->integrinStrength();

			double temp_integrins = pC->custom_data["integrin"];

			temp_a *= dens * temp_integrins;
			
			tmp_r -= temp_a;
		}
		
		/////////////////////////////////////////////////////////////////
		if(tmp_r==0)
			return;
		tmp_r/=distance;

		pC->velocity += tmp_r * pC->displacement;
	}
}

void custom_update_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	
	//First check the neighbors in my current voxel
	for( auto neighbor : pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()] )
	{
		pCell->add_potentials( neighbor );
	}

	pCell->state.simple_pressure = 0.0; 
	pCell->state.neighbors.clear(); // new 1.8.0

	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		pCell->add_potentials(*neighbor);
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			pCell->add_potentials(*neighbor);
		}
	}

	int ecm_index = BioFVM::microenvironment.find_density_index("ECM_collagen");
	if ( ecm_index >= 0 ){
		add_ecm_interaction( pCell, ecm_index, pCell->get_current_mechanics_voxel_index() );
	}
/*
	if (pCell->custom_data["freezed"] > 2){
		return ;
	}
*/
	
	pCell->update_motility_vector(dt);
	//std::cout << pCell->state.simple_pressure << " \n ";
	pCell->velocity += phenotype.motility.motility_vector;
	
	return; 
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ 
	
	
	
	return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

bool read_custom_microenvironment_from_matlab( std::string mat_filename )
{
	std::cout << std::endl << "Attempting to load the microenvironment from " << mat_filename << " ... " << std::endl; 

	std::vector< std::vector<double> > mat = read_matlab( mat_filename ); 

	// row 0 : x
	// row 1 : y
	// row 2 : z 
	// row 3 : vol 
	// row 4-n : substrate 
	int num_rows = mat.size(); 
	int num_cols = mat[0].size(); 

	int number_of_mat_voxels = num_cols; 
	int number_of_mat_substrates = num_rows - 3 -1; 

	if( number_of_mat_substrates != microenvironment.number_of_densities() )
	{
		std::cout << "Error reading microenvironment from " << mat_filename << "! ";  
		std::cout << "Expected " << microenvironment.number_of_densities() << " substrates but only detected "
			<< number_of_mat_substrates << std::endl; 
		return false; 
	}

	if( number_of_mat_voxels != microenvironment.number_of_voxels() )
	{
		std::cout << "Error reading microenvironment from " << mat_filename << "! ";  
		std::cout << "Expected " << microenvironment.number_of_voxels() << " voxels but only detected "
			<< number_of_mat_voxels << std::endl; 
		return false; 
	}

	std::cout << num_rows << std::endl;

	std::cout << num_cols << std::endl;
	int ecm_index = microenvironment.find_density_index("ECM_collagen");
	for( int n=0 ; n < number_of_mat_voxels ; n++ )
	{
		std::cout << microenvironment(n) << std::endl;

		std::cout << mat[5][n] << std::endl;
		for( int k=4; k < num_rows ; k++ )
		{ 
			microenvironment(n)[k-4] = mat[k][n]; 
			}
		//microenvironment.density_vector(n)[ecm_index] = mat[5][n];
	}

	std::cout << "done!" << std::endl << std::endl; 
	return true; 
}

void add_low_cut(){

	int ecm_index = microenvironment.find_density_index("ECM_collagen");
	double low_cut = PhysiCell::parameters.doubles("low_cut");
	for( int n=0 ; n < microenvironment.number_of_voxels() ; n++ )
	{
		if(microenvironment.density_vector(n)[ecm_index] < low_cut){

			microenvironment.density_vector(n)[ecm_index] = 0;
		}
		
	}
}