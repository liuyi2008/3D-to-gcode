#ifndef main_h
#define main_h
#include <iostream>  
#include <OpenMesh/Core/IO/MeshIO.hh>  
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>  
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh\Core\Utils\Property.hh> 
#include "GL\glut.h"

#include "GL\freeglut.h"
#include <math.h>
#include <Windows.h>
#include <string>
#include <vector>
#include <algorithm>
//#include "smooth_algo.hh"

using namespace std;
struct MyTraits : public OpenMesh::DefaultTraits
{
	HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
	VertexTraits
	{
		int some_additional_index;
	};
	FaceTraits{
		int cd_add_index;
	};
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

void initGL();

void myReshape(GLint w, GLint h);

#endif
