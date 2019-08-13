#ifndef slice_h
#define slice_h
//#include <windows.h>
//#include <stdlib.h> 
#include <iostream>  
#include <OpenMesh/Core/IO/MeshIO.hh>  
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>  
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh\Core\Utils\Property.hh> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include<GL/glut.h>
#include<stdio.h>
#include "GL\glut.h"
#include "GL\freeglut.h"
using namespace std;

struct point {

	double x, y;

};

extern vector<vector<point> >model;


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

//vector�ṹ�����ڴ�������뽻��
class culine {
public:
	double StLine[3] = { 0.0 };
	double EdLine[3] = { 0.0 };
	double CrossPoint[3] = { 0.0 };
};
//���ڴ洢���潻��ͽ�������

void readfile(string file);// ��ȡ�ļ��ĺ���

void BoundingBox();//ͨ���򵥵ı�������ģ�Ͱ�Χ�У����ȡ�����޹أ�

bool IntersectPlane(MyMesh::Point pt, MyMesh::Point pnorm, MyMesh::Point *pilist, int &pnum);//���溯��

void findIntersect();

#endif

