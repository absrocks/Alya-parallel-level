
#ifdef ALYA_GMSH
#include <gmshc.h>
#endif
#include <stdio.h>
#include <stdlib.h>


void addModel(char *name_model,
		    int  tag){
  // Some variables declarations
  int     ierr;

#ifdef ALYA_GMSH
  // --> 1. Add a gmsh model and their respective entities
  gmshModelAdd(name_model, &ierr);  
  // Add entities in all three dimensions 
  gmshModelAddDiscreteEntity(1, tag, NULL, 0, &ierr); // ndimension = 1
  gmshModelAddDiscreteEntity(2, tag, NULL, 0, &ierr); // ndimension = 2
  gmshModelAddDiscreteEntity(3, tag, NULL, 0, &ierr); // ndimension = 3

#endif
}
  
/*  
------------------------------------------------------
  addMesh subroutine
------------------------------------------------------
  ELEMENT TYPE INFORMATION 
  
    type element	 1 is actually "Line 2" 
    type element	 2 is actually "Triangle 3" 
    type element	 3 is actually "Quadrilateral 4" 
    type element	 4 is actually "Tetrahedron 4" 
    type element	 5 is actually "Hexahedron 8" 
    type element	 6 is actually "Prism 6" 
    type element	 7 is actually "Pyramid 5" 
    type element	 8 is actually "Line 3" 
    type element	 9 is actually "Triangle 6" 
    type element        10 is actually "Quadrilateral 9" 
  ...
------------------------------------------------------
  INPUT VARIABLES

  Nodal information
  1. nnodes:                     numer of nodes,  
  2. node_tags[nnodes]:          array of node tags,             
  3. node_coordinates[3*nnodes]: array of nodal coordinates. 

  Element information
  4. ntypes:                                                           number of element types
  5. nelements_bytype[ntypes]:                                         array with the number of elements by type:  
  6. element_tags[nelements_bytype[1] +
                  nelements_bytype[2] +  ...
		  nelements_bytype[ntypes]:                            arrays of element tags order by type:
  7. elements_nodes[nelements_bytype[1]*nnodes_byelement[1] +
                    nelements_bytype[2]*nnodes_byelement[1] + ...
                    nelements_bytype[ntypes]*nnodes_byelement[ntypes]: arrays of elemental nodes order by type,
     where nnodes_byelement[itype] is the nodes of element by type that is NOT given as an input is deduced.
------------------------------------------------------
  EXAMPLE FOR ADDING A MESH

  // Define nodes tags and coordinates
  size_t nodes_tags[9]          = {1,  2,  3,  4,  5,   6,  7,  8, 9};  
  double nodes_coordinates[9*3] = {0., 0., 0.,
				   1., 0., 0.,
				   2., 0., 0.,
				   0., 1., 0.,
				   1., 1., 0.,
				   2., 1., 0.,
				   0., 2., 0.,
				   1., 2., 0.,
				   2., 2., 0.};
  // Define types: quadrilateral and  triangles
  size_t types_tags[2]       = {3, 2};  
  // Number of elements for each type
  size_t nelements_bytype[2] = {2, 4}; 
  // All elements tags order by type
  size_t elements_tags[6]    = {1, 2, 3, 4, 5, 6};
  // All elements nodes definitions order by type
  size_t elements_nodes[20]  = {1, 2, 5, 4,
			        2, 3, 6, 5,
				4, 8, 7,
				4, 5, 8,
				5, 9, 8,
				5, 6, 9};
  // Add input mesh information to gmsh
  addMesh(9, // nnodes = 9
	  nodes_tags,	 
	  nodes_coordinates,
	  2, // ntypes = 2
	  types_tags,	       
	  nelements_bytype, 
	  elements_tags,	 
	  elements_nodes);
*/
void addMesh(size_t   nnodes,
	     size_t  *nodes_tags,	 
	     double  *nodes_coordinates,
	     size_t   ntypes,	     
	     size_t  *types_tags, 	       
	     size_t  *nelements_bytype, 
	     size_t  *elements_tags,	 
	     size_t  *elements_nodes,

	     char    *name_model,  //input parameter
	     int      ndime,       //input/output parameter  if "ndime" < 0 then add all types of elments ans save the maximum dimension in "ndime"
	     int      tag){        //input parameter         if "ndime" > 0 then add only the elments wirtha dimension equals to "ndime"
                                     
  // Some variables declarations
  int     ierr;
  int     ielement;
  int     ielement_nodes;
  int     dim;
  int     order;
  int     nnodes_byelement;
  double *nodeCoord;  
  size_t  nodeCoord_n;  
  int     numPrimaryNodes;
  char   *elementName;
  size_t  itype;
  size_t  i;
#ifdef ALYA_GMSH
  
  // --> 1. Add nodes to gmsh
  //printf(">> ADDING A MESH  \n");
  gmshModelMeshAddNodes(3,
			tag,
			nodes_tags,
			nnodes,
			nodes_coordinates,
			nnodes*3,
			NULL,
			0,
			&ierr);
  // printf(">>   Error adding %zu nodes = %i  \n", nnodes, ierr);

  ielement       = 0;  // index of elements
  ielement_nodes = 0;  // index of elemental nodes
  // Loop over element types  
  for(itype = 0; itype < ntypes; itype++) {
    // Obtain element type information
    gmshModelMeshGetElementProperties(types_tags[itype], // type of element
				      &elementName,
				      &dim,
				      &order,					
				      &nnodes_byelement, // number of nodes by element for a type
				      &nodeCoord,
				      &nodeCoord_n,
				      &numPrimaryNodes,
				      &ierr);
    
    if (ndime < 0 || dim == ndime) { // Save only the elements equal to "ndime", if ndime < 0 save all of them
      // --> 2. Add elements to gmsh by types
      gmshModelMeshAddElementsByType(tag,
				     types_tags[itype],
				     &elements_tags[ielement],
				     nelements_bytype[itype],
				     &elements_nodes[ielement_nodes],
				     nelements_bytype[itype]*nnodes_byelement,
				     &ierr);
      //printf(">>   Error adding %zu elements of type \"%s\" = %i  \n", nelements_bytype[itype], elementName, ierr);
    }    
    ielement       = ielement       + nelements_bytype[itype];      
    ielement_nodes = ielement_nodes + nelements_bytype[itype]*nnodes_byelement;

    free(nodeCoord);
    free(elementName);    
  }

#endif
  
}



void addViewModel(size_t  body_nelements,
		  size_t *body_elements_tags,
		  double *mesh_size_input,
		  char   *model_name,
		  int    *sf_view){

  // Some variables declarations
  int            ierr;
  size_t        *ones_array;
  double       **mesh_size;
  const double **mesh_size_const;
  int i;

#ifdef ALYA_GMSH
  // --> 1. Add a new view 
  *sf_view = gmshViewAdd("mesh size field",-1,&ierr);
  // printf(">>   Error adding view = %i  \n", ierr);

  // --> 2.1  Allocate space for size field
  ones_array = (size_t *) malloc(body_nelements * sizeof(size_t));
  mesh_size  = (double **)malloc(body_nelements * sizeof(double*));
  for(i = 0; i < body_nelements; i++){
    mesh_size[i] = (double *)malloc(1 * sizeof(double));
  }
  // --> 2.2 Save data in size field with the input data
  for(i = 0; i < body_nelements; i++){
    mesh_size [i][0] = mesh_size_input[i];
    ones_array[i]    = 1;
  }
  mesh_size_const = (const double **)mesh_size;
  
  // --> 3. Add model data to the plrevious view.
  //        Here, we pass the input size mesh field to gmsh
  gmshViewAddModelData(*sf_view,
		       0,
		       model_name,
		       "ElementData",
		       body_elements_tags,
		       body_nelements,
		       mesh_size_const,
		       ones_array,		       
		       body_nelements,
		       0.0, -1, 0, &ierr );
  
  // printf(">>   Error adding model = %i  \n", ierr);  
  gmshViewWrite(*sf_view,"sf_view.pos",0,&ierr);

  // --> 3.1  Free pointer for store the space for size field
  for(i = 0; i < body_nelements; i++){
    free(mesh_size[i]);
  }
  free(mesh_size);
  free(ones_array);
#endif
}          


/*  
------------------------------------------------------
  adaptMesh subroutine
------------------------------------------------------
*/
void adaptMesh(int sf_view,
	       int ndime){
  // Some variables declarations
  int ierr;

  //printf(">> ADAPTING MESH  \n");  
  // --> 1. Create a surface geometry in gmsh from a boundary mesh previously added
  //        Also create an empty body geomtry from the surface geometry previously added  
  int curves_tags[1] = {2};

#ifdef ALYA_GMSH 
  
  if (ndime == 2){
    
    int boundary_tag     = gmshModelGeoAddCurveLoop(curves_tags,1,-1,&ierr); // connect 1D boundary mesh with the gmsh geometry
    int boundary_tags[1] = {boundary_tag};
    int body_tag         = gmshModelGeoAddPlaneSurface(boundary_tags,1,-1,&ierr);

    gmshOptionSetNumber("Mesh.Algorithm", 5, &ierr); // 2d
    
  } else if (ndime == 3){
     
    int boundary_tag     = gmshModelGeoAddSurfaceLoop(curves_tags,1,-1,&ierr); // connect 2D boundary mesh with the gmsh geometry
    int boundary_tags[1] = {boundary_tag};    
    int body_tag         = gmshModelGeoAddVolume(boundary_tags,1,-1,&ierr);

    gmshOptionSetNumber("Mesh.Algorithm3D", 10, &ierr); // 3d
    
  }
  
  gmshModelGeoSynchronize(&ierr);
  //printf(">>   Error synchronizing = %i  \n", ierr);

  // Add a background field
  
  int bg_field = gmshModelMeshFieldAdd("PostView",-1,&ierr);
  
  gmshModelMeshFieldSetNumber(bg_field, "ViewTag", sf_view,&ierr);
  
  gmshModelMeshFieldSetAsBackgroundMesh(bg_field,&ierr);
  
  // Generate the new adapted mesh
  gmshModelMeshGenerate(ndime,&ierr);
  
  //printf(">>   Error generating the adapted mesh in %i-D = %i  \n", ndime,ierr);
#endif
  
}          




/*  
------------------------------------------------------
  getMesh subroutine
------------------------------------------------------
*/
void getMesh(size_t  *new_nnodes,           
	     size_t **new_nodes_tags,       
	     double **new_nodes_coordinates,
	     size_t  *new_ntypes,           
	     size_t **new_types_tags,       
	     size_t **new_nelements_bytype, 
	     size_t **new_elements_tags,    
	     size_t **new_elements_nodes,
	     int     ndime,
	     int     tag,
	     size_t *new_nelements_total)
{
  // Some variables declarations
  size_t  coord_n, parametricCoord_n;
  double *parametricCoord;
  int     ierr;  

  size_t **elementTags;	   
  size_t   elementTags_nn; 
  size_t **nodeTags;	   
  size_t  *nodeTags_n;	   
  size_t   nodeTags_nn;
  int     *new_types_tags_tem;

  size_t *new_types_tags_tem2;
  size_t nelements;
  size_t nelements_nodes;
  size_t *elements_tags_tem;
  size_t *elements_nodes_tem;
  size_t itype;
  size_t ielement;
  size_t inode;

#ifdef ALYA_GMSH
 
  // printf(">> GETTING A MESH  \n");    
  // --> 1. Get all the mesh nodes from gmsh with a dimension equals to ndime
  gmshModelMeshGetNodes(new_nodes_tags, new_nnodes,
			new_nodes_coordinates, &coord_n,
			&parametricCoord, &parametricCoord_n,
			-1,
			-1,
			0,
			0,
			&ierr);
  // printf(">>   Error getting %zu nodes = %i  \n", *new_nnodes, ierr);
  
  // --> 2. Get all the mesh elements from gmsh with a dimension equals to ndime
  gmshModelMeshGetElements(&new_types_tags_tem, new_ntypes,
			   &elementTags, new_nelements_bytype, &elementTags_nn,
			   &nodeTags, &nodeTags_n, &nodeTags_nn,
			   ndime,
			   tag,
			   &ierr);
  // for(size_t itype = 0; itype < *new_ntypes; itype++) printf(">>   Error getting %zu elements  = %i  \n", *(*new_nelements_bytype+itype), ierr);
  
  // Transform types tags from int to size_t
  new_types_tags_tem2 = calloc(*new_ntypes, sizeof(size_t));
  for(itype = 0; itype < *new_ntypes; itype++) {
    new_types_tags_tem2[itype]  = (size_t)new_types_tags_tem[itype];
  }  

  // --> 3. Determine the total number of elements tags and nodes considering all types
  *new_nelements_total = 0;
  for(itype = 0; itype < *new_ntypes; itype++) {      
    *new_nelements_total = *new_nelements_total + *(*new_nelements_bytype+itype);
  }  
  nelements_nodes = 0;
  for(itype = 0; itype < *new_ntypes; itype++) {      
    nelements_nodes = nelements_nodes + nodeTags_n[itype];    
  }

  // --> 4. Store elemetal tags in temporal array
  elements_tags_tem = calloc(*new_nelements_total, sizeof(size_t));
  nelements= 0;
  // Loop over types
  for(itype = 0; itype < *new_ntypes; itype++) {
    // Loop over elemental nodes by type
    for(ielement = 0; ielement < *(*new_nelements_bytype+itype); ielement++) {
      // Save the elemwntal tags
      elements_tags_tem[itype*nelements + ielement]  = elementTags[itype][ielement];
    }
    nelements = nelements + *(*new_nelements_bytype+itype);
  }
  
  // --> 5. Store elemetal nodes in a temporal array
  elements_nodes_tem = calloc(nelements_nodes, sizeof(size_t));  
  nelements_nodes = 0;
  // Loop over types  
  for(itype = 0; itype < *new_ntypes; itype++) {
    // Loop over elemental nodes by type      
    for(inode = 0; inode < nodeTags_n[itype]; inode++) {
      // Save the elemtal nodes      
      elements_nodes_tem[itype*nelements_nodes + inode]  = nodeTags[itype][inode];  
    }
    nelements_nodes = nelements_nodes + nodeTags_n[itype];
  }
  
  // --> 6. Store elemental tags and nodes from the temporal arrays to the input ones
  *new_elements_nodes = elements_nodes_tem;
  *new_elements_tags  = elements_tags_tem;
  *new_types_tags     = new_types_tags_tem2;

  free(parametricCoord);
  for(itype = 0; itype<*new_ntypes; itype++){
     free(elementTags[itype]);
     free(nodeTags[itype]);
  }
  free(elementTags);	
  free(nodeTags); 
  free(nodeTags_n); 
  free(new_types_tags_tem);

#endif
} 



 
void freeArrays(size_t **new_nodes_tags,
		double **new_nodes_coordinates,
		size_t **new_types_tags,
		size_t **new_nelements_bytype,
		size_t **new_elements_tags,	 
		size_t **new_elements_nodes){
  free(*new_nodes_tags);
  free(*new_nodes_coordinates);
  free(*new_types_tags);
  free(*new_nelements_bytype);
  free(*new_elements_tags);	 
  free(*new_elements_nodes);  
}


void reMeshing(int      ndime,
	       
	       size_t   nnodes,
	       size_t  *nodes_tags,	 
	       double  *nodes_coordinates,
	       size_t   ntypes,	     
	       size_t  *types_tags,	       
	       size_t  *nelements_bytype, 
	       size_t  *elements_tags,	 
	       size_t  *elements_nodes,

	       double   min_size,
	       double   max_size,	       
	       double  *mesh_size_input,

	       size_t  *new_nnodes,
	       size_t **new_nodes_tags,	 
	       double **new_nodes_coordinates,
	       size_t  *new_ntypes,	     
	       size_t **new_types_tags,	       
	       size_t **new_nelements_bytype, 
	       size_t **new_elements_tags,	 
	       size_t **new_elements_nodes){
  
  // Variables declarations
  int     tag;
  char   *model_name;
    
  int     ierr;
  char   *argv;

  size_t  body_nnodes;		   
  size_t *body_nodes_tags;	    
  double *body_nodes_coordinates;
  size_t  body_ntypes;		  
  size_t *body_types_tags;	  
  size_t *body_nelements_bytype; 
  size_t *body_elements_tags;	  
  size_t *body_elements_nodes;	      
  size_t  body_nelements;
  
  int     sf_view;

    
#ifdef ALYA_GMSH
  
  //***********************************************************************  
  // --> 1. Initialize gmsh
  //***********************************************************************
  gmshInitialize(0, &argv, 1, &ierr);
  tag = 1; model_name = "input";  
  addModel(model_name,tag);
  
  // Do not use boundary sizes to fix minimum size
  gmshOptionSetNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0, &ierr);
  // Size field maximun and minimum 
  gmshOptionSetNumber("Mesh.CharacteristicLengthMin", min_size, &ierr);
  gmshOptionSetNumber("Mesh.CharacteristicLengthMax", max_size, &ierr);
  
  //***********************************************************************
  // --> 2. Add input mesh information to gmsh
  //***********************************************************************

  addMesh(nnodes, 
	  nodes_tags,	 
	  nodes_coordinates,
	  ntypes, 
	  types_tags,	       
	  nelements_bytype, 
	  elements_tags,	 
	  elements_nodes,
	  model_name,
          -1, // ndime = -1 (add all types of elements)
	  tag);

  // Write input mesh
  // gmshWrite("input.msh", &ierr);

  //***********************************************************************
  // --> 3. Remesh the domain according to input size field using gmsh
  //***********************************************************************  
  // --> 3.1 Get input element information from gmsh without include the boundary
  getMesh(&body_nnodes,
  	  &body_nodes_tags,
  	  &body_nodes_coordinates,
  	  &body_ntypes,
  	  &body_types_tags,
  	  &body_nelements_bytype,
  	  &body_elements_tags,
  	  &body_elements_nodes,
  	  ndime, 
	  tag,
  	  &body_nelements);
    
  // --> 3.2 Save mesh size field data adding a view model data in gmsh and create a new global model
  addViewModel(body_nelements,
	       body_elements_tags,
	       mesh_size_input,
	       model_name,
	       &sf_view);  
  // --> 3.3 Create a new model for the adapated mesh
  tag = 2; model_name = "output";  
  addModel(model_name,tag);

  // --> 3.3 Add only input boundary elements to the new model
  //         Problem (1): ALL input nodes will be included
  addMesh(nnodes, 
  	  nodes_tags,
  	  nodes_coordinates,
  	  ntypes,
  	  types_tags,
  	  nelements_bytype,
  	  elements_tags,
  	  elements_nodes,
  	  model_name,
          ndime-1,  // add only the boundary elements
	  tag);
  
  // --> 3.4 Adapat mesh with the information given in mesh_size
  //     Problem (2): input BOUNDARY nodes are used in the new adapted mesh
  adaptMesh(sf_view,ndime);
  
  // --> 3.5 Preserve only the nodes related with some element in the adapted mesh
  //         Solution to (1-2): keep only the input nodes related with the boundary
  //         Problem(3):        nodes will be not sorted in a sequential way.  
  gmshModelMeshReclassifyNodes(&ierr);
  // Renumerate nodes sequently
  // Solution to (3): sort nodes list sequently    
  gmshModelMeshRenumberNodes(&ierr);

  // Writhe adapated mesh
  // gmshWrite("adapted.msh",&ierr);
  
  //***********************************************************************	         
  // --> 4. Get output mesh information from gmsh (without boundary mesh)
  //***********************************************************************
  //     Problem (4): elements will be not sorted in a sequential way
  size_t new_nelements;
  getMesh(new_nnodes,	       
	  new_nodes_tags,	  
	  new_nodes_coordinates,
	  new_ntypes,		  
	  new_types_tags,	  
	  new_nelements_bytype, 
	  new_elements_tags,	  
	  new_elements_nodes,
	  ndime, // only body mesh
	  tag+1, // new tag value for the adapated mesh (+1 of the last tag) 
	  &new_nelements);
  //***********************************************************************


  //***********************************************************************	         
  // --> 5. Obtain an elemental sequance list 
  //***********************************************************************  
  // --> 5.1 Create a new modesl 
  tag = 10; model_name = "final";  
  addModel(model_name,tag);
  // --> 5.2 Add the adapted mesh to gmsh without including boundary mesh    
  addMesh(*new_nnodes,
  	  *new_nodes_tags,
  	  *new_nodes_coordinates,
  	  *new_ntypes,
  	  *new_types_tags,
  	  *new_nelements_bytype,
  	  *new_elements_tags,
  	  *new_elements_nodes,
  	  "final",
          ndime,   // add only the body elements)
  	  tag);

  // Renumerate eleemnts from 1 and sequently
  // Solution to (4): sort elemental list sequently    
  gmshModelMeshRenumberElements(&ierr);

  // --> 5.3 Free memory 
  freeArrays(new_nodes_tags,
	     new_nodes_coordinates,
	     new_types_tags,
	     new_nelements_bytype,
	     new_elements_tags,	 
	     new_elements_nodes);
  // --> 5.4 Get the adapated mesh arrays without boundary mesh and sequently sorted 
  getMesh(new_nnodes,
  	  new_nodes_tags,
  	  new_nodes_coordinates,
  	  new_ntypes,
  	  new_types_tags,
  	  new_nelements_bytype,
  	  new_elements_tags,
  	  new_elements_nodes,
  	  ndime, // only body mesh
  	  tag, // new tag value for the adapated mesh (+1 of the last tag)
  	  &new_nelements);
  
  // Write the adapated mesh without boundary mesh 
  gmshWrite("output.msh",&ierr);

  // Finalize and free some pointers     
  gmshFinalize(&ierr);
  free(body_nodes_tags);	    
  free(body_nodes_coordinates);
  free(body_types_tags);	  
  free(body_nelements_bytype); 
  free(body_elements_tags);	  
  free(body_elements_nodes);

#else

  // If ALYA_GMSH is not actived, copy the input mesh as the output one   
  size_t  inode;
  double *tem_nodes_coordinates;  
  tem_nodes_coordinates = (double *) malloc(3*nnodes * sizeof(double));
  for(inode = 0; inode < 3*nnodes; inode++) {
    tem_nodes_coordinates[inode] = nodes_coordinates[inode];
  }
  
  *new_nnodes            = nnodes;     
  *new_nodes_tags        = nodes_tags;
  *new_nodes_coordinates = tem_nodes_coordinates;  
  *new_ntypes          	 = ntypes;          
  *new_types_tags        = types_tags;      
  *new_nelements_bytype  = nelements_bytype;
  *new_elements_tags     = elements_tags;   
  *new_elements_nodes	 = elements_nodes;
  
#endif
  
}

/* *********************************************************************** */
/* printf(" nnodes : %zu \n",*new_nnodes); */
/* printf(" ***************************** \n"); */
/* for(size_t i = 0; i < *new_nnodes; i++) { */
/*   printf(" node tag  %zu: %zu \n",i,*(*new_nodes_tags+i)); */
/*   printf(" new_nodes %zu: %f %f %f \n",i,*(*new_nodes_coordinates+3*i),*(*new_nodes_coordinates+3*i+1),*(*new_nodes_coordinates+3*i+2)); */
/* } */
/* printf(" ***************************** \n"); */
/* printf(" ntypes : %zu \n",*new_ntypes ); */
/* for(size_t i = 0; i < *new_ntypes; i++) { */
/*   printf(" nelements_bytype %zu: %zu \n",i,*(*new_nelements_bytype+i)); */

/* } */
/* printf(" ***************************** \n"); */
/* for(int i = 0; i < 3; i++) { */
/*   printf(" element tag   %d: %zu \n",i,*(*new_elements_tags+i) ); */
/*   printf(" element nodes %d: %zu %zu \n",i,*(*new_elements_nodes+2*i),*(*new_elements_nodes+2*i+1)); */
/* } */
/* for(int i = 0; i < 8; i++) { */
/*   printf(" element tag   %d: %zu \n",i,*(*new_elements_tags+i) ); */
/*   printf(" element nodes %d: %zu %zu %zu \n",i+3,*(*new_elements_nodes+3*i+6),*(*new_elements_nodes+3*i+1+6),*(*new_elements_nodes+3*i+2+6)); */
/* } */
/* printf(" ***************************** \n"); */
/* *********************************************************************** */
