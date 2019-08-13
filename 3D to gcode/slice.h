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

//vector结构体用于储存截面与交点
class culine {
public:
	double StLine[3] = { 0.0 };
	double EdLine[3] = { 0.0 };
	double CrossPoint[3] = { 0.0 };
};
//用于存储截面交点和交点数量

void readfile(string file);// 读取文件的函数

void BoundingBox();//通过简单的遍历点获得模型包围盒（与获取截面无关）

bool IntersectPlane(MyMesh::Point pt, MyMesh::Point pnorm, MyMesh::Point *pilist, int &pnum);//截面函数

void findIntersect();

#endif

