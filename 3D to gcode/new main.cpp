#include <iostream>  
#include <OpenMesh/Core/IO/MeshIO.hh>  
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>  
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh\Core\Utils\Property.hh> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include<stdio.h>
#include <gl/glut.h>  
//#include<GL/glut.h>
//#include "GL\glut.h"
//#include "GL\freeglut.h"
using namespace std;

struct point {

	double x, y;

}p,p1;

typedef unsigned char boolean;

#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)

vector<point>coord;
vector<point>polypoint; //一层的轮廓线坐标
vector<point>tripoint;  //一层的变化线 线段坐标 成对
vector<point>new_tripoint;  //新的一层的变化线 线段坐标 成对
vector<point>interpoint; //交点坐标
vector<vector<point> >model; //slice  整个模型的，分层的轮廓线坐标，一维层数，二维此层的轮廓线点坐标
vector<vector<point> >modelfill(model); //完成的model填充线，一维层数，二维此层的变化线 线段坐标

const string file_1 = "test.obj";  
const string file_2 = "buddha_head.obj";
const string file_3 = "wawa.obj";
const string file_4 = "opener.obj";
const string file_5 = "test2.obj";
const string file_6 = "bunnyr.obj";
const string file_7 = "half buddha.obj";
const string file_8 = "Mouse.obj";
const string file_9 = "Mouse.stl";
int currentfile = 1;

const float z = 0.20;   //层厚
const float t = 0.060;   //挤出率

//与调节视角有关
double xs = 0;
double ys = 0;
double zs = 0;

//与调节视点有关
double sx = 0;
double sy = 0;
double sz = 0;

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
const float ERR = 0.001;

//文件读取有关的
MyMesh mesh;

GLuint showFaceList, showWireList, showCutList;
int showstate = 1;
bool showFace = true;
bool showWire = false;
bool showCut = true;
bool showFlatlines = false;

//用于包围盒与截面
double Bmax_x, Bmax_y, Bmax_z, Bmin_x, Bmin_y, Bmin_z, px, py, pz;
const int INF_formistake = 0x3f3f3f3f;
int mid_n = 0; int mid_m = 0;
int countt = 0;

MyMesh::Point plist[3000][1000];
int numofcut[3000];

// 读取文件的函数
void readfile(string file) {
	// 请求顶点法线 vertex normals
	mesh.request_vertex_normals();
	//如果不存在顶点法线，则报错 
	if (!mesh.has_vertex_normals())
	{
		cout << "错误：标准定点属性 “法线”不存在" << endl;
		return;
	}
	// 如果有顶点发现则读取文件 
	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, file, opt))
	{
		cout << "无法读取文件:" << file << endl;
		return;
	}
	else cout << "成功读取文件:" << file << endl;
	cout << endl; // 为了ui显示好看一些
				  //如果不存在顶点法线，则计算出
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		// 通过面法线计算顶点法线
		mesh.request_face_normals();
		// mesh计算出顶点法线
		mesh.update_normals();
		// 释放面法线
		mesh.release_face_normals();
	}
}

//通过简单的遍历点获得模型包围盒（与获取截面无关）
void BoundingBox() {
	MyMesh::Point pt;
	int st = 0;
	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
		MyMesh::VertexHandle vh_i = *it;
		pt = mesh.point(vh_i);
		px = pt.data()[0];
		py = pt.data()[1];
		pz = pt.data()[2];
		if (st == 0) {
			Bmax_x = Bmin_x = px;
			Bmax_y = Bmin_y = py;
			Bmax_z = Bmin_z = pz;
			st++;
		}
		else {
			if (px > Bmax_x)Bmax_x = px; else if (px < Bmin_x)Bmin_x = px;
			if (py > Bmax_y)Bmax_y = py; else if (py < Bmin_y)Bmin_y = py;
			if (pz > Bmax_z)Bmax_z = pz; else if (pz < Bmin_z)Bmin_z = pz;
		}
	}
	//	Bmax_x += 0.2; Bmax_y += 0.2; Bmax_z += 0.2;
	//	Bmin_x -= 0.2; Bmin_y -= 0.2; Bmin_z -= 0.2;
	//printf("%f %f %f %f %f %f\n", Bmax_x, Bmax_y, Bmax_z, Bmin_x, Bmin_y, Bmin_z);
}

//截面函数
bool IntersectPlane(MyMesh::Point pt, MyMesh::Point pnorm, MyMesh::Point *pilist, int &pnum) {
	//参数包括 pt，pnorm，*pilist，pnum[]  具体函数原理 见 截面算法.docx
	int starte, ne, ne1, nf, num;
	MyMesh::Point vt1, vt2;
	//MyMesh::Face f1;
	MyMesh::HalfedgeHandle nhe;
	MyMesh::FaceHandle nhf;
	float d1, d2, sd1, sd2;
	bool *flag, suc;
	float dist, mind = 1.0e+8;

	pnum = 0;
	sd1 = sd2 = -10000;
	int esize = mesh.n_halfedges();
	flag = new bool[esize];

	suc = false;


	for (MyMesh::HalfedgeIter it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) {

		MyMesh::HalfedgeHandle hh = *it;
		int id = hh.idx();
		flag[id] = false;

		auto fromVertex = mesh.from_vertex_handle(hh);
		auto toVertex = mesh.to_vertex_handle(hh);
		vt1 = mesh.point(fromVertex);
		vt2 = mesh.point(toVertex);
		//printf("$ %.3f %.3f $\n", vt1.data()[0],vt2.data()[0]);
		d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
		d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);

		if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0) || d1 > 0 && d2 == 0 || d1 == 0 && d2 > 0))
		{
			flag[id] = true;

			vt1.data()[0] = vt1.data()[0] - pt.data()[0];
			vt1.data()[1] = vt1.data()[1] - pt.data()[1];
			vt1.data()[2] = vt1.data()[2] - pt.data()[2];       // point date minus point date 
			dist = vt1.data()[0] * vt1.data()[0] + vt1.data()[1] * vt1.data()[1] + vt1.data()[2] * vt1.data()[2];
			if (dist < mind) {
				nhe = hh;
				mind = dist;
				ne = id;                    //  printf("ne:  %d  \n", ne);
				suc = true;
			}
		}
	}

	if (!suc) {
		delete[]flag;
		return false;
	}

	starte = ne;
	num = 0;

	suc = false;

	nhf = mesh.face_handle(nhe);

	while (!suc) {
		//printf("%%%%");	

		auto fromVertex = mesh.from_vertex_handle(nhe);
		auto toVertex = mesh.to_vertex_handle(nhe);

		vt1 = mesh.point(fromVertex);
		vt2 = mesh.point(toVertex);

		d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
		d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);
		//printf("$$$%lf %lf \n", d1, d2);
		if ((sd1 == d1) && (sd2 == d2)) {
			flag[ne] = false;
		}
		sd1 = d1; sd2 = d2;
		pilist[num].data()[0] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0];
		pilist[num].data()[1] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1];
		pilist[num].data()[2] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[2] - vt1.data()[2]) + vt1.data()[2];
		//printf("$$$%lf %lf %lf %lf %lf %lf\n", vt1.data()[0], vt1.data()[1], vt1.data()[2], vt2.data()[0], vt2.data()[1], vt2.data()[2]);
		num++;



		int nn = 0;
		do {
			for (auto it = mesh.fh_begin(nhf); it != mesh.fh_end(nhf); ++it) {

				MyMesh::HalfedgeHandle halfnow = *it;
				nn++;
				if (nn > 3) {
					starte = ne;
					flag[ne] = false;
					break;
				}
				const int ne1 = halfnow.idx();

				if (flag[ne1] == false || ne == ne1) continue;
				nn = 0;
				MyMesh::VertexHandle fromV = mesh.from_vertex_handle(halfnow);
				MyMesh::VertexHandle toV = mesh.to_vertex_handle(halfnow);

				vt1 = mesh.point(fromV);
				vt2 = mesh.point(toV);

				d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
					+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
				d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
					+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);

				if (((d1 >= 0 && d2 <= 0) || (d1 <= 0 && d2 >= 0)) && fabs(d1 - d2) > ERR) //ERR的值未知 
				{

					MyMesh::HalfedgeHandle halfnext = mesh.opposite_halfedge_handle(halfnow);

					nhf = mesh.face_handle(halfnext);


					int ne2 = halfnext.idx();

					flag[ne1] = flag[ne2] = false;


					if (nhf.idx() == -1)
					{
						starte = ne;
						flag[ne] = false;
						break;
					}
					ne = ne2;
					// compute the intersecting point		
					pilist[num].data()[0] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0];
					pilist[num].data()[1] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1];
					pilist[num].data()[2] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[2] - vt1.data()[2]) + vt1.data()[2];
					//printf("##%lf %lf %lf\n", pilist[num].data()[0], pilist[num].data()[1], pilist[num].data()[2]);
					num++;

					break;
				}
			}
		} while (ne != starte);

		suc = true;

		for (auto it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) {

			MyMesh::HalfedgeHandle hh = *it;

			int id = hh.idx();

			if (flag[id] == true) {

				ne = id;

				starte = ne;
				nhe = hh;
				nhf = mesh.face_handle(nhe);
				if (nhf.idx() == -1)
				{
					flag[ne] = false;

					continue;
				}
				pilist[num].data()[0] = -10000;
				pilist[num].data()[1] = -10000;
				pilist[num].data()[2] = -10000;
				num++;


				suc = false;
				break;
			}
		}


	};

	pnum = num;
	//printf("%d\n", pnum);

	delete[]flag;

	return true;
}


void findIntersect() {
	countt = 0;
	for (double i = Bmin_z; i < Bmax_z; ) {
		MyMesh::Point pt;
		MyMesh::Point pilist[1000];
		MyMesh::Normal vf(0, 0, 1);
		int pnum = 0;
		float Xport = 0;
		pt.data()[0] = 0; pt.data()[1] = 0; pt.data()[2] = i;
		if (IntersectPlane(pt, vf, pilist, pnum)) {
			numofcut[countt] = pnum;
			if (Bmin_x * Bmax_x < 0) {

				if (i < 0) {
					Xport = 100 * (abs(i - Bmin_x) / (abs(Bmax_x) + abs(Bmin_x)));
				}
				else {
					Xport = 100 * ((i + abs(Bmin_x)) / (abs(Bmax_x) + abs(Bmin_x)));
				}
			}
			else {
				Xport = 100 * (abs(i - Bmin_x) / abs(abs(Bmax_x) - abs(Bmin_x)));
			}
			//cout << "######"<< Xport << endl;
			//printf("-------------------------------------------------\n");
			if (Xport < 18) {
				//printf("第%d个截面环\t    %.4f%%    头部\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport / 10);
			}
			else if (Xport < 25) {
				//printf("第%d个截面环\t    %.4f%%    颈部\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport / 10);
			}
			else if (Xport < 35) {
				//printf("第%d个截面环\t    %.4f%%    前胸\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport / 10);
			}
			else if (Xport < 77) {
				//printf("第%d个截面环\t    %.4f%%    躯干\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport / 10);
			}
			else if (Xport < 96) {
				//printf("第%d个截面环\t    %.4f%%    臀部\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport / 10);
			}
			else {
				//printf("第%d个截面环\t    %.4f%%    尾部\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport / 10);
			}
			int j; int numdiff = 0; MyMesh::Point diff = pilist[0];
			for (j = 0; j < pnum; j++) {
				plist[countt][j] = pilist[j];
				if (j == pnum - 1 || ((pilist[j + 1].data()[0] == -10000) && (pilist[j + 3].data()[0] != -10000))) {

					if ((pilist[j].data()[2] - diff.data()[2]) < 0.001) {
						numdiff++;
						continue;
					}//cout << diff.data()[2] << pilist[j].data()[2] << endl;
				}
				if (pilist[j].data()[0] == -10000) {
					if (j < pnum - 1)diff = pilist[j + 1];
					//printf("-- %d th:(   EX  , EX  , EX   ) \t\t", j + 1 - numdiff);
				}
				else
				{
					//printf("-- %d th:( %.4f, %.4f, %.4f )   \t", j + 1 - numdiff, pilist[j].data()[0], pilist[j].data()[1], pilist[j].data()[2]);
					p = { pilist[j].data()[0] ,pilist[j].data()[1] };
					coord.push_back(p);
				}
				if ((j + numdiff) % 2) {
					//printf("\n");
				}
			}

			/*
			int j;
			for (j = 0; j < pnum; j++) {
				plist[countt][j] = pilist[j];
				if (pilist[j].data()[0] == -10000) {
					printf("-- %d th:(   EX  , EX  , EX   ) \t\t",j+1);
				}else printf("-- %d th:( %.4f, %.4f, %.4f )   \t", j+1, pilist[j].data()[0], pilist[j].data()[1], pilist[j].data()[2]);
				if (j % 2) {
					printf("\n");
				}
			}
			*/
			if (j % 2) {
				//printf("\n");
			}
			//printf("-------------------------------------------------\n");
			model.push_back(coord);
			coord.clear();
			countt++;
		}
		/*
		//选择性的对关键部位提高截面密度
		if ((i > Bmin_x + abs(Bmax_x - Bmin_x)*0.75 && i < Bmin_x + abs(Bmax_x - Bmin_x)*0.87) ||
			(i<Bmin_x + abs(Bmax_x - Bmin_x)*0.45 && i>Bmin_x + abs(Bmax_x - Bmin_x)*0.33) ||
			(i > Bmin_x + abs(Bmax_x - Bmin_x)*0.92)) {
			i += 20;
		}
		else {
			i += 5;
		}
		*/
		i += z;
	}
}

//***********************************************************************以上切片一下填充****************************************************************
//*************判断线段与线段的关系函数*************************
double xmult(point p1, point p2, point p0) {
	return (p1.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p1.y - p0.y);
}

//判点是否在线段上,包括端点
int dot_online_in(point p, point l1, point l2) {
	return zero(xmult(p, l1, l2)) && (l1.x - p.x)*(l2.x - p.x) < eps && (l1.y - p.y)*(l2.y - p.y) < eps;
}

//判两点在线段同侧,点在线段上返回0
int same_side(point p1, point p2, point l1, point l2) {
	return xmult(l1, p1, l2)*xmult(l1, p2, l2) > eps;
}

//判两直线平行
int parallel(point u1, point u2, point v1, point v2) {
	return zero((u1.x - u2.x)*(v1.y - v2.y) - (v1.x - v2.x)*(u1.y - u2.y));
}

//判三点共线
int dots_inline(point p1, point p2, point p3) {
	return zero(xmult(p1, p2, p3));
}

//判两线段相交,包括端点和部分重合
int intersect_in(point u1, point u2, point v1, point v2) {
	if (!dots_inline(u1, u2, v1) || !dots_inline(u1, u2, v2))
		return !same_side(u1, u2, v1, v2) && !same_side(v1, v2, u1, u2);
	return dot_online_in(u1, v1, v2) || dot_online_in(u2, v1, v2) || dot_online_in(v1, u1, u2) || dot_online_in(v2, u1, u2);
}

//计算两线段交点,请判线段是否相交(同时还是要判断是否平行!)
point intersection(point u1, point u2, point v1, point v2) {
	point ret = u1;
	double t = ((u1.x - v1.x)*(v1.y - v2.y) - (u1.y - v1.y)*(v1.x - v2.x))
		/ ((u1.x - u2.x)*(v1.y - v2.y) - (u1.y - u2.y)*(v1.x - v2.x));
	ret.x += (u2.x - u1.x)*t;
	ret.y += (u2.y - u1.y)*t;
	return ret;
}

//主函数，有交点存交点，没交点无操作
void inter(point u1, point u2, point v1, point v2)  
{
	point ans;
	if (parallel(u1, u2, v1, v2) || !intersect_in(u1, u2, v1, v2))
	{

	}
	else {
		ans = intersection(u1, u2, v1, v2);
		interpoint.push_back(ans);
		//printf("交点为:(%lf,%lf)", ans.x, ans.y);
	}

}

//*************判断线段与线段的关系函数*************************

int InOrOutPolygon(point a)  //判断点是否在多边形内
{
	double x0 = a.x;
	double y0 = a.y;
	int crossings = 0;
	int n = polypoint.size();
	for (int i = 0; i < n; i++)
	{
		// 点在两个x之间 且以点垂直y轴向上做射线
		double slope = (polypoint[(i + 1 + n) % n].y - polypoint[i].y) / (polypoint[(i + 1 + n) % n].x - polypoint[i].x);
		boolean cond1 = (polypoint[i].x <= x0) && (x0 < polypoint[(i + 1 + n) % n].x);
		boolean cond2 = (polypoint[(i + 1 + n) % n].x <= x0) && (x0 < polypoint[i].x);
		boolean above = (y0 < slope * (x0 - polypoint[i].x) + polypoint[i].y);
		if ((cond1 || cond2) && above) crossings++;
	}
	return (crossings % 2 != 0);    //返 回 值:  0:外 1:内
}

point trans(point p, float b, float c)  //坐标平移  这里似乎不应该缩放，只需要平移
{
	//a 缩放倍数,b 左右平移,c 上下平移
	p = { p.x + b , p.y + c };
	return p;
}

//point trans2(point p)//调整轮廓线
//{
//	int a = 1;  //缩放
//	int b = -250;  //左右平移
//	int c = 0;  //上下平移
//	p = { a * p.x + b , a * p.y + c };
//	return p;
//}

point mindistance(point a, point b, point c) //在b 和 c中返回离a 最近的点
{
	int dis1 = pow((a.x - b.x), 2) + pow((a.y - b.y), 2);
	int dis2 = pow((a.x - c.x), 2) + pow((a.y - c.y), 2);
	if (dis1 < dis2)return b;
	else return c;
}

double distance(point a, point b)
{
	return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2));
}

void triangle(int i)  //第i层，从1开始，由i和w一起决定点A、B、C的位置  1层
{
	float a = 30.0f;// 三角形边长 这里来控制缩放

	float width = 0.4; //喷头的宽度,或者变化的速率，现在看来就是每次偏移的距离，和真实碰头宽度配合决定了角度，即变化率

	int m = 0;  //行
	int n = 0;  //列
	//************************处理层数i,使之映射在区间0到4*sqrt(3)*0.333*a里*****************************
	while (i - 1 > 4 * sqrt(3)*0.333*a / width)
	{
		i = i - int(4 * sqrt(3)*0.333*a / width);
	}
	int condition;
	if (i - 1 < sqrt(3)*0.333*a / width)
	{
		condition = 1;
	}
	else if (i - 1 < 2 * sqrt(3)*0.333*a / width)
	{
		condition = 2;
	}
	else if (i - 1 < 3 * sqrt(3)*0.333*a / width)
	{
		condition = 3;
	}
	else if (i - 1 < 4 * sqrt(3)*0.333*a / width)
	{
		condition = 4;
	}
	//************************处理层数i,使之映射在区间0到4*sqrt(3)*0.333*a里*****************************
	switch (condition)
	{
	case 1:
	{
		for (m; m < (120 * 2) / (sqrt(3)*a); m++)  //阶段一：w*i 在 0——sqrt(3)/3 *a
		{
			for (n; n < 170 / a; n++)
			{
				point A = { m*0.5*a + n * a + 0.5*a ,                         m*0.5*sqrt(3)*a + sqrt(3) * 0.5 * a - width * (i - 1) };

				point B = { m*0.5*a + n * a + 0.5*sqrt(3)*width*(i - 1) ,       m*0.5*sqrt(3)*a + 0.5*width * (i - 1) };

				point C = { m*0.5*a + n * a + a - 0.5*sqrt(3)*width*(i - 1) ,   m*0.5*sqrt(3)*a + 0.5*width * (i - 1) };

				point D = { m*0.5*a + n * a + 0.5*a ,                         m*0.5*sqrt(3)*a + 0.5*sqrt(3)*a };

				point E = { m * 0.5 * a + n * a + 0 ,                         m*0.5*sqrt(3)*a + 0 };

				point F = { m*0.5*a + n * a + a ,                             m*0.5*sqrt(3)*a + 0 };

				tripoint.push_back(trans(E,-50,0)); tripoint.push_back(trans(B, -50, 0));
				tripoint.push_back(trans(B, -50, 0)); tripoint.push_back(trans(C, -50, 0));
				tripoint.push_back(trans(C, -50, 0)); tripoint.push_back(trans(F, -50, 0));
				//tripoint.push_back(trans(F,-50,0)); tripoint.push_back(trans(C,-50,0));
				tripoint.push_back(trans(C, -50, 0)); tripoint.push_back(trans(A, -50, 0));
				tripoint.push_back(trans(A, -50, 0)); tripoint.push_back(trans(D, -50, 0));
				//tripoint.push_back(trans(D,-50,0)); tripoint.push_back(trans(A,-50,0));
				tripoint.push_back(trans(A, -50, 0)); tripoint.push_back(trans(B, -50, 0));
				//tripoint.push_back(trans(B,-50,0)); tripoint.push_back(trans(E,-50,0));

			}
			n = 0;
		}
		break;
	}
	case 2:
	{
		for (m; m < (120 * 2) / (sqrt(3)*a); m++)   //阶段二：w*i在 sqrt(3)/3 *a——2 * sqrt(3)/3 *a
		{
			for (n; n < 170 / a; n++)
			{
				point A = { m*0.5*a + n * a + 0.5*a ,                           m*0.5*sqrt(3)*a + width * (i - 1) - sqrt(3)*0.333*a + sqrt(3)*0.167*a };

				point B = { m*0.5*a + n * a + a - sqrt(3)*0.5*width*(i - 1) ,     m*0.5*sqrt(3)*a + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };

				point C = { m*0.5*a + n * a + sqrt(3)*0.5*width*(i - 1) ,       m*0.5*sqrt(3)*a + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };

				point D = { m*0.5*a + n * a + 0.5*a ,                           m*0.5*sqrt(3)*a + 0.5*sqrt(3)*a };

				point E = { m * 0.5 * a + n * a + 0 ,                           m*0.5*sqrt(3)*a + 0 };

				point F = { m*0.5*a + n * a + a ,                               m*0.5*sqrt(3)*a + 0 };

				tripoint.push_back(trans(E, -50, 0)); tripoint.push_back(trans(B, -50, 0));
				tripoint.push_back(trans(B, -50, 0)); tripoint.push_back(trans(C, -50, 0));
				tripoint.push_back(trans(C, -50, 0)); tripoint.push_back(trans(F, -50, 0));
				//tripoint.push_back(trans(F,-50,0)); tripoint.push_back(trans(C,-50,0));
				tripoint.push_back(trans(C, -50, 0)); tripoint.push_back(trans(A, -50, 0));
				tripoint.push_back(trans(A, -50, 0)); tripoint.push_back(trans(D, -50, 0));
				//tripoint.push_back(trans(D,-50,0)); tripoint.push_back(trans(A,-50,0));
				tripoint.push_back(trans(A, -50, 0)); tripoint.push_back(trans(B, -50, 0));
				//tripoint.push_back(trans(B,-50,0)); tripoint.push_back(trans(E,-50,0));

			}
			n = 0;
		}
		break;
	}
	case 3:
	{
		i = i - 2 * sqrt(3)*0.333*a / width;
		for (m; m < (120 * 2) / (sqrt(3)*a); m++)  //阶段三：和阶段一的区别只是Y坐标全变成负
		{
			for (n; n < 170 / a; n++)
			{
				point A = { m*0.5*a + n * a + 0.5*a ,                         m*0.5*sqrt(3)*a - sqrt(3) * 0.5 * a + width * (i - 1) };

				point B = { m*0.5*a + n * a + 0.5*sqrt(3)*width*(i - 1) ,       m*0.5*sqrt(3)*a - 0.5*width * (i - 1) };

				point C = { m*0.5*a + n * a + a - 0.5*sqrt(3)*width*(i - 1) ,   m*0.5*sqrt(3)*a - 0.5*width * (i - 1) };

				point D = { m*0.5*a + n * a + 0.5*a ,                          m*0.5*sqrt(3)*a - 0.5*sqrt(3)*a };

				point E = { m * 0.5 * a + n * a + 0 ,                          m*0.5*sqrt(3)*a + 0 };

				point F = { m*0.5*a + n * a + a ,                              m*0.5*sqrt(3)*a + 0 };

				tripoint.push_back(trans(E, -50, 0)); tripoint.push_back(trans(B, -50, 0));
				tripoint.push_back(trans(B, -50, 0)); tripoint.push_back(trans(C, -50, 0));
				tripoint.push_back(trans(C, -50, 0)); tripoint.push_back(trans(F, -50, 0));
				//tripoint.push_back(trans(F,-50,0)); tripoint.push_back(trans(C,-50,0));
				tripoint.push_back(trans(C, -50, 0)); tripoint.push_back(trans(A, -50, 0));
				tripoint.push_back(trans(A, -50, 0)); tripoint.push_back(trans(D, -50, 0));
				//tripoint.push_back(trans(D,-50,0)); tripoint.push_back(trans(A,-50,0));
				tripoint.push_back(trans(A, -50, 0)); tripoint.push_back(trans(B, -50, 0));
				//tripoint.push_back(trans(B,-50,0)); tripoint.push_back(trans(E,-50,0));

			}
			n = 0;
		}
		break;
	}
	case 4:
	{
		i = i - 2 * sqrt(3)*0.333*a / width;
		for (m; m < (120 * 2) / (sqrt(3)*a); m++)  //阶段四：和阶段二一样
		{
			for (n; n < 170 / a; n++)
			{
				point A = { m*0.5*a + n * a + 0.5*a ,                           m*0.5*sqrt(3)*a - width * (i - 1) + sqrt(3)*0.333*a - sqrt(3)*0.167*a };

				point B = { m*0.5*a + n * a + a - sqrt(3)*0.5*width*(i - 1) ,     m*0.5*sqrt(3)*a - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };

				point C = { m*0.5*a + n * a + sqrt(3)*0.5*width*(i - 1) ,       m*0.5*sqrt(3)*a - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };

				point D = { m*0.5*a + n * a + 0.5*a ,                           m*0.5*sqrt(3)*a - 0.5*sqrt(3)*a };

				point E = { m * 0.5 * a + n * a + 0 ,                           m*0.5*sqrt(3)*a + 0 };

				point F = { m*0.5*a + n * a + a ,                               m*0.5*sqrt(3)*a + 0 };

				tripoint.push_back(trans(E, -50, 0)); tripoint.push_back(trans(B, -50, 0));
				tripoint.push_back(trans(B, -50, 0)); tripoint.push_back(trans(C, -50, 0));
				tripoint.push_back(trans(C, -50, 0)); tripoint.push_back(trans(F, -50, 0));
				//tripoint.push_back(trans(F,-50,0)); tripoint.push_back(trans(C,-50,0));
				tripoint.push_back(trans(C, -50, 0)); tripoint.push_back(trans(A, -50, 0));
				tripoint.push_back(trans(A, -50, 0)); tripoint.push_back(trans(D, -50, 0));
				//tripoint.push_back(trans(D,-50,0)); tripoint.push_back(trans(A,-50,0));
				tripoint.push_back(trans(A, -50, 0)); tripoint.push_back(trans(B, -50, 0));
				//tripoint.push_back(trans(B,-50,0)); tripoint.push_back(trans(E,-50,0));

			}
			n = 0;
		}
		break;
	}

	default: printf("有问题");
		break;
	}

}

//void triangle(int i)  //第i层，从1开始，由i和w一起决定点A、B、C的位置         两层
//{
//	float a = 20.0f;// 三角形边长 这里来控制缩放
//
//	float width = 0.4; //喷头的宽度,或者变化的速率，现在看来就是每次偏移的距离，和真实碰头宽度配合决定了角度，即变化率
//
//	float x = a + 2 * sqrt(3)*0.33*width;  //x坐标右移单位
//	float y = 0.5*sqrt(3)*a + width;  //y坐标上移单位
//
//	int m = 0;  //行
//	int n = 0;  //列
//	//************************处理层数i,使之映射在区间0到4*sqrt(3)*0.333*a里*****************************
//	while (i - 1 > 4 * sqrt(3)*0.333*a / width)
//	{
//		i = i - int(4 * sqrt(3)*0.333*a / width);
//	}
//	int condition;
//	if (i - 1 < sqrt(3)*0.333*a / width)
//	{
//		condition = 1;
//	}
//	else if (i - 1 < 2 * sqrt(3)*0.333*a / width)
//	{
//		condition = 2;
//	}
//	else if (i - 1 < 3 * sqrt(3)*0.333*a / width)
//	{
//		condition = 3;
//	}
//	else if (i - 1 < 4 * sqrt(3)*0.333*a / width)
//	{
//		condition = 4;
//	}
//	//************************处理层数i,使之映射在区间0到4*sqrt(3)*0.333*a里*****************************
//	switch (condition)
//	{
//	case 1:
//	{
//		// for (m; m < 1; m++)  //阶段一：w*i 在 0——sqrt(3)/3 *a
//		// {
//		   //for (n; n < 1; n++)
//		   //{
//		for (m; m < (120 * 2) / (sqrt(3)*a); m++)  //阶段一：w*i 在 0——sqrt(3)/3 *a
//		{
//			for (n; n < 170 / a; n++)
//			{
//				//point A = { m*0.5*a + n * a + 0.5*a ,                         m*0.5*sqrt(3)*a + sqrt(3) * 0.5 * a - width * (i - 1) };
//
//				//point B = { m*0.5*a + n * a + 0.5*sqrt(3)*width*(i - 1) ,       m*0.5*sqrt(3)*a + 0.5*width * (i - 1) };
//
//				//point C = { m*0.5*a + n * a + a - 0.5*sqrt(3)*width*(i - 1) ,   m*0.5*sqrt(3)*a + 0.5*width * (i - 1) };
//
//				//point D = { m*0.5*a + n * a + 0.5*a ,                         m*0.5*sqrt(3)*a + 0.5*sqrt(3)*a };
//
//				//point E = { m * 0.5 * a + n * a + 0 ,                         m*0.5*sqrt(3)*a + 0 };
//
//				//point F = { m*0.5*a + n * a + a ,                             m*0.5*sqrt(3)*a + 0 };
//
//				point A = { m*0.5*x + n * x + 0.5*a ,                         m * y + sqrt(3) * 0.5 * a - width * (i - 1) };
//
//				point A1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*width ,    m*y + sqrt(3) * 0.5 * a - width * (i - 1) + 0.5*width };
//
//				point A2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*width,      m*y + sqrt(3) * 0.5 * a - width * (i - 1) + 0.5*width };
//
//				point B = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) ,       m*y + 0.5*width * (i - 1) };
//
//				point B1 = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) - 0.33*sqrt(3)*width ,       m*y + 0.5*width * (i - 1) };
//
//				point B2 = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) - 0.167*sqrt(3)*width ,       m*y + 0.5*width * (i - 1) - 0.5*width };
//
//				point C = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) ,   m*y + 0.5*width * (i - 1) };
//
//				point C1 = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) + 0.167*sqrt(3)*width ,   m*y + 0.5*width * (i - 1) - 0.5*width };
//
//				point C2 = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) + 0.33*sqrt(3)*width ,   m*y + 0.5*width * (i - 1) };
//
//				point D1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*width,                         m*y + 0.5*sqrt(3)*a + 0.5*width };
//
//				point D2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*width,                         m*y + 0.5*sqrt(3)*a + 0.5*width };
//
//				point E1 = { m * 0.5 * x + n * x + 0 - 0.33*sqrt(3)*width ,                         m*y + 0 };
//
//				point E2 = { m * 0.5 * x + n * x + 0 - 0.167*sqrt(3)*width ,                         m*y + 0 - 0.5*width };
//
//				point F1 = { m*0.5*x + n * x + a + 0.33*sqrt(3)*width,                             m*y + 0 };
//
//				point F2 = { m*0.5*x + n * x + a + 0.167*sqrt(3)*width,                             m*y + 0 - 0.5*width };
//
//				tripoint.push_back(trans(B1)); tripoint.push_back(trans(A1));
//				tripoint.push_back(trans(A1)); tripoint.push_back(trans(D1));
//				tripoint.push_back(trans(D1)); tripoint.push_back(trans(D2));
//				tripoint.push_back(trans(D2)); tripoint.push_back(trans(A2));
//				tripoint.push_back(trans(A2)); tripoint.push_back(trans(C2));
//				tripoint.push_back(trans(C2)); tripoint.push_back(trans(F1));
//				tripoint.push_back(trans(F1)); tripoint.push_back(trans(F2));
//				tripoint.push_back(trans(F2)); tripoint.push_back(trans(C1));
//				tripoint.push_back(trans(C1)); tripoint.push_back(trans(B2));
//				tripoint.push_back(trans(B2)); tripoint.push_back(trans(E2));
//				tripoint.push_back(trans(E2)); tripoint.push_back(trans(E1));
//				tripoint.push_back(trans(E1)); tripoint.push_back(trans(B1));
//				//tripoint.push_back(trans(B1)); tripoint.push_back(trans(B));
//
//				tripoint.push_back(trans(B)); tripoint.push_back(trans(A));
//				tripoint.push_back(trans(A)); tripoint.push_back(trans(C));
//				tripoint.push_back(trans(C)); tripoint.push_back(trans(B));
//
//			}
//			n = 0;
//		}
//		break;
//	}
//	case 2:
//	{
//		for (m; m < (120 * 2) / (sqrt(3)*a); m++)   //阶段二：w*i在 sqrt(3)/3 *a——2 * sqrt(3)/3 *a
//		{
//			for (n; n < 170 / a; n++)
//			{
//				point A = { m*0.5*x + n * x + 0.5*a ,                           m*y + width * (i - 1) - sqrt(3)*0.333*a + sqrt(3)*0.167*a };
//
//				point A1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*width,     m*y + width * (i - 1) - sqrt(3)*0.333*a + sqrt(3)*0.167*a + 0.5*width };
//
//				point A2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*width,     m*y + width * (i - 1) - sqrt(3)*0.333*a + sqrt(3)*0.167*a + 0.5*width };
//
//				point B = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) ,     m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };
//
//				point B1 = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) - 0.33*sqrt(3)*width ,     m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };
//
//				point B2 = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) - 0.167*sqrt(3)*width ,     m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) - 0.5*width };
//
//				point C = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) ,       m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };
//
//				point C1 = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) + 0.167*sqrt(3)*width ,       m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) - 0.5*width };
//
//				point C2 = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) + 0.33*sqrt(3)*width,       m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };
//
//				point D1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*width,                         m*y + 0.5*sqrt(3)*a + 0.5*width };
//
//				point D2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*width,                         m*y + 0.5*sqrt(3)*a + 0.5*width };
//
//				point E1 = { m * 0.5 * x + n * x + 0 - 0.33*sqrt(3)*width ,                         m*y + 0 };
//
//				point E2 = { m * 0.5 * x + n * x + 0 - 0.167*sqrt(3)*width ,                         m*y + 0 - 0.5*width };
//
//				point F1 = { m*0.5*x + n * x + a + 0.33*sqrt(3)*width,                             m*y + 0 };
//
//				point F2 = { m*0.5*x + n * x + a + 0.167*sqrt(3)*width,                             m*y + 0 - 0.5*width };
//
//				tripoint.push_back(trans(B1)); tripoint.push_back(trans(A1));
//				tripoint.push_back(trans(A1)); tripoint.push_back(trans(D1));
//				tripoint.push_back(trans(D1)); tripoint.push_back(trans(D2));
//				tripoint.push_back(trans(D2)); tripoint.push_back(trans(A2));
//				tripoint.push_back(trans(A2)); tripoint.push_back(trans(C2));
//				tripoint.push_back(trans(C2)); tripoint.push_back(trans(F1));
//				tripoint.push_back(trans(F1)); tripoint.push_back(trans(F2));
//				tripoint.push_back(trans(F2)); tripoint.push_back(trans(C1));
//				tripoint.push_back(trans(C1)); tripoint.push_back(trans(B2));
//				tripoint.push_back(trans(B2)); tripoint.push_back(trans(E2));
//				tripoint.push_back(trans(E2)); tripoint.push_back(trans(E1));
//				tripoint.push_back(trans(E1)); tripoint.push_back(trans(B1));
//				//tripoint.push_back(trans(B1)); tripoint.push_back(trans(B));
//
//				tripoint.push_back(trans(B)); tripoint.push_back(trans(A));
//				tripoint.push_back(trans(A)); tripoint.push_back(trans(C));
//				tripoint.push_back(trans(C)); tripoint.push_back(trans(B));
//
//			}
//			n = 0;
//		}
//		break;
//	}
//	case 3:
//	{
//		i = i - 2 * sqrt(3)*0.333*a / width;
//		for (m; m < (120 * 2) / (sqrt(3)*a); m++)  //阶段三：和阶段一的区别只是Y坐标全变成负
//		{
//			for (n; n < 170 / a; n++)
//			{
//				point A = { m*0.5*x + n * x + 0.5*a ,                         m * y - sqrt(3) * 0.5 * a + width * (i - 1) };
//
//				point A1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*width ,    m*y - sqrt(3) * 0.5 * a + width * (i - 1) - 0.5*width };
//
//				point A2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*width,      m*y - sqrt(3) * 0.5 * a + width * (i - 1) - 0.5*width };
//
//				point B = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) ,       m*y - 0.5*width * (i - 1) };
//
//				point B1 = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) - 0.33*sqrt(3)*width ,       m*y - 0.5*width * (i - 1) };
//
//				point B2 = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) - 0.167*sqrt(3)*width ,       m*y - 0.5*width * (i - 1) + 0.5*width };
//
//				point C = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) ,   m*y - 0.5*width * (i - 1) };
//
//				point C1 = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) + 0.167*sqrt(3)*width ,   m*y - 0.5*width * (i - 1) + 0.5*width };
//
//				point C2 = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) + 0.33*sqrt(3)*width ,   m*y - 0.5*width * (i - 1) };
//
//				point D1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*width,                         m*y - 0.5*sqrt(3)*a - 0.5*width };
//
//				point D2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*width,                         m*y - 0.5*sqrt(3)*a - 0.5*width };
//
//				point E1 = { m * 0.5 * x + n * x + 0 - 0.33*sqrt(3)*width ,                         m*y + 0 };
//
//				point E2 = { m * 0.5 * x + n * x + 0 - 0.167*sqrt(3)*width ,                         m*y + 0 + 0.5*width };
//
//				point F1 = { m*0.5*x + n * x + a + 0.33*sqrt(3)*width,                             m*y + 0 };
//
//				point F2 = { m*0.5*x + n * x + a + 0.167*sqrt(3)*width,                             m*y + 0 + 0.5*width };
//
//				tripoint.push_back(trans(B1)); tripoint.push_back(trans(A1));
//				tripoint.push_back(trans(A1)); tripoint.push_back(trans(D1));
//				tripoint.push_back(trans(D1)); tripoint.push_back(trans(D2));
//				tripoint.push_back(trans(D2)); tripoint.push_back(trans(A2));
//				tripoint.push_back(trans(A2)); tripoint.push_back(trans(C2));
//				tripoint.push_back(trans(C2)); tripoint.push_back(trans(F1));
//				tripoint.push_back(trans(F1)); tripoint.push_back(trans(F2));
//				tripoint.push_back(trans(F2)); tripoint.push_back(trans(C1));
//				tripoint.push_back(trans(C1)); tripoint.push_back(trans(B2));
//				tripoint.push_back(trans(B2)); tripoint.push_back(trans(E2));
//				tripoint.push_back(trans(E2)); tripoint.push_back(trans(E1));
//				tripoint.push_back(trans(E1)); tripoint.push_back(trans(B1));
//				//tripoint.push_back(trans(B1)); tripoint.push_back(trans(B));
//
//				tripoint.push_back(trans(B)); tripoint.push_back(trans(A));
//				tripoint.push_back(trans(A)); tripoint.push_back(trans(C));
//				tripoint.push_back(trans(C)); tripoint.push_back(trans(B));
//
//			}
//			n = 0;
//		}
//		break;
//	}
//	case 4:
//	{
//		i = i - 2 * sqrt(3)*0.333*a / width;
//		for (m; m < (120 * 2) / (sqrt(3)*a); m++)  //阶段四：和阶段二一样
//		{
//			for (n; n < 170 / a; n++)
//			{
//
//				point A = { m*0.5*x + n * x + 0.5*a ,                           m*y - width * (i - 1) + sqrt(3)*0.333*a - sqrt(3)*0.167*a };
//
//				point A1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*width,     m*y - width * (i - 1) + sqrt(3)*0.333*a - sqrt(3)*0.167*a - 0.5*width };
//
//				point A2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*width,     m*y - width * (i - 1) + sqrt(3)*0.333*a - sqrt(3)*0.167*a - 0.5*width };
//
//				point B = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) ,     m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };
//
//				point B1 = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) - 0.33*sqrt(3)*width ,     m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };
//
//				point B2 = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) - 0.167*sqrt(3)*width ,     m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) + 0.5*width };
//
//				point C = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) ,       m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };
//
//				point C1 = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) + 0.167*sqrt(3)*width ,       m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) + 0.5*width };
//
//				point C2 = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) + 0.33*sqrt(3)*width,       m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };
//
//				point D1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*width,                         m*y - 0.5*sqrt(3)*a - 0.5*width };
//
//				point D2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*width,                         m*y - 0.5*sqrt(3)*a - 0.5*width };
//
//				point E1 = { m * 0.5 * x + n * x + 0 - 0.33*sqrt(3)*width ,                         m*y + 0 };
//
//				point E2 = { m * 0.5 * x + n * x + 0 - 0.167*sqrt(3)*width ,                         m*y + 0 + 0.5*width };
//
//				point F1 = { m*0.5*x + n * x + a + 0.33*sqrt(3)*width,                             m*y + 0 };
//
//				point F2 = { m*0.5*x + n * x + a + 0.167*sqrt(3)*width,                             m*y + 0 + 0.5*width };
//
//				tripoint.push_back(trans(B1)); tripoint.push_back(trans(A1));
//				tripoint.push_back(trans(A1)); tripoint.push_back(trans(D1));
//				tripoint.push_back(trans(D1)); tripoint.push_back(trans(D2));
//				tripoint.push_back(trans(D2)); tripoint.push_back(trans(A2));
//				tripoint.push_back(trans(A2)); tripoint.push_back(trans(C2));
//				tripoint.push_back(trans(C2)); tripoint.push_back(trans(F1));
//				tripoint.push_back(trans(F1)); tripoint.push_back(trans(F2));
//				tripoint.push_back(trans(F2)); tripoint.push_back(trans(C1));
//				tripoint.push_back(trans(C1)); tripoint.push_back(trans(B2));
//				tripoint.push_back(trans(B2)); tripoint.push_back(trans(E2));
//				tripoint.push_back(trans(E2)); tripoint.push_back(trans(E1));
//				tripoint.push_back(trans(E1)); tripoint.push_back(trans(B1));
//				//tripoint.push_back(trans(B1)); tripoint.push_back(trans(B));
//
//				tripoint.push_back(trans(B)); tripoint.push_back(trans(A));
//				tripoint.push_back(trans(A)); tripoint.push_back(trans(C));
//				tripoint.push_back(trans(C)); tripoint.push_back(trans(B));
//
//			}
//			n = 0;
//		}
//		break;
//	}
//
//	default: printf("有问题");
//		break;
//	}
//
//
//
//}

//void triangle(int i)  //第i层，从1开始，由i和w一起决定点A、B、C的位置         三层
//{
//	float a = 30.0f;// 三角形边长 这里来控制缩放 当三角形厚度增加，不可不考虑壁厚，这里是内三角形边长
//
//	float width = 0.4; //上层和下层的偏移的距离，和真实喷头宽度配合决定了角度，即变化率，范围(0,一个喷宽)
//	//这个变量控制线与线之间的距离，在下面的句意中所有路径间距0.5*width，在实际打印中喷头宽度要大于0.5*width,打双层时要等于0.5*width，打3层时等于0.25*width
//
//	float offset = 1.6;  //同层不同路径之间的距离，等于一个喷头宽度
//
//	float x = a + 2 * sqrt(3)*0.33*offset + 0.25*offset;  //x坐标右移单位 这里是外三角形边长
//	float y = 0.5*sqrt(3)*a + offset + 0.25*offset;  //y坐标上移单位
//
//	int m = 0;  //行
//	int n = 0;  //列
//	//************************处理层数i,使之映射在区间0到4*sqrt(3)*0.333*a里*****************************
//	while (i - 1 > 4 * sqrt(3)*0.333*a / width)
//	{
//		i = i - int(4 * sqrt(3)*0.333*a / width);
//	}
//	int condition;
//	if (i - 1 < sqrt(3)*0.333*a / width)
//	{
//		condition = 1;
//	}
//	else if (i - 1 < 2 * sqrt(3)*0.333*a / width)
//	{
//		condition = 2;
//	}
//	else if (i - 1 < 3 * sqrt(3)*0.333*a / width)
//	{
//		condition = 3;
//	}
//	else if (i - 1 < 4 * sqrt(3)*0.333*a / width)
//	{
//		condition = 4;
//	}
//	//************************处理层数i,使之映射在区间0到4*sqrt(3)*0.333*a里*****************************
//	switch (condition)
//	{
//	case 1:
//	{
//		// for (m; m < 1; m++)  //阶段一：w*i 在 0——sqrt(3)/3 *a
//		// {
//		   //for (n; n < 1; n++)
//		   //{
//		for (m; m < (120 * 2) / (sqrt(3)*a); m++)  //阶段一：w*i 在 0——sqrt(3)/3 *a
//		{
//			for (n; n < 170 / a; n++)
//			{
//				point A3 = { m*0.5*x + n * x + 0.5*a ,                         m*y + sqrt(3) * 0.5 * a - width * (i - 1) + 0.5*offset };
//
//				point B3 = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) - 0.25*sqrt(3)*offset ,       m*y + 0.5*width * (i - 1) - 0.25*offset };
//
//				point C3 = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) + 0.25*sqrt(3)*offset,   m*y + 0.5*width * (i - 1) - 0.25*offset };
//
//				point D3 = { m*0.5*x + n * x + 0.5*a ,                         m*y + 0.5*sqrt(3)*a + 0.5*offset };
//
//				point E3 = { m * 0.5 * x + n * x + 0 - 0.25*sqrt(3)*offset,                         m*y + 0 - 0.25*offset };
//
//				point F3 = { m*0.5*x + n * x + a + 0.25*sqrt(3)*offset,                             m*y + 0 - 0.25*offset };
//
//				point A = { m*0.5*x + n * x + 0.5*a ,                         m * y + sqrt(3) * 0.5 * a - width * (i - 1) };
//
//				point A1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*offset ,    m*y + sqrt(3) * 0.5 * a - width * (i - 1) + 0.5*offset };
//
//				point A2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*offset,      m*y + sqrt(3) * 0.5 * a - width * (i - 1) + 0.5*offset };
//
//				point B = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) ,       m*y + 0.5*width * (i - 1) };
//
//				point B1 = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) - 0.33*sqrt(3)*offset ,       m*y + 0.5*width * (i - 1) };
//
//				point B2 = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) - 0.167*sqrt(3)*offset ,       m*y + 0.5*width * (i - 1) - 0.5*offset };
//
//				point C = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) ,   m*y + 0.5*width * (i - 1) };
//
//				point C1 = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) + 0.167*sqrt(3)*offset ,   m*y + 0.5*width * (i - 1) - 0.5*offset };
//
//				point C2 = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) + 0.33*sqrt(3)*offset ,   m*y + 0.5*width * (i - 1) };
//
//				point D1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*offset,                         m*y + 0.5*sqrt(3)*a + 0.5*offset };
//
//				point D2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*offset,                         m*y + 0.5*sqrt(3)*a + 0.5*offset };
//
//				point E1 = { m * 0.5 * x + n * x + 0 - 0.33*sqrt(3)*offset ,                         m*y + 0 };
//
//				point E2 = { m * 0.5 * x + n * x + 0 - 0.167*sqrt(3)*offset ,                         m*y + 0 - 0.5*offset };
//
//				point F1 = { m*0.5*x + n * x + a + 0.33*sqrt(3)*offset,                             m*y + 0 };
//
//				point F2 = { m*0.5*x + n * x + a + 0.167*sqrt(3)*offset,                             m*y + 0 - 0.5*offset };
//
//
//				tripoint.push_back(trans(E3, -50, 0)); tripoint.push_back(trans(E1, -50, 0));
//				tripoint.push_back(trans(E1, -50, 0)); tripoint.push_back(trans(B1, -50, 0));
//				tripoint.push_back(trans(B1, -50, 0)); tripoint.push_back(trans(A1, -50, 0));
//				tripoint.push_back(trans(A1, -50, 0)); tripoint.push_back(trans(D1, -50, 0));
//				tripoint.push_back(trans(D1, -50, 0)); tripoint.push_back(trans(D2, -50, 0));
//				tripoint.push_back(trans(D2, -50, 0)); tripoint.push_back(trans(A2, -50, 0));
//				tripoint.push_back(trans(A2, -50, 0)); tripoint.push_back(trans(C2, -50, 0));
//				tripoint.push_back(trans(C2, -50, 0)); tripoint.push_back(trans(F1, -50, 0));
//				tripoint.push_back(trans(F1, -50, 0)); tripoint.push_back(trans(F2, -50, 0));
//				tripoint.push_back(trans(F2, -50, 0)); tripoint.push_back(trans(C1, -50, 0));
//				tripoint.push_back(trans(C1, -50, 0)); tripoint.push_back(trans(B2, -50, 0));
//				tripoint.push_back(trans(B2, -50, 0)); tripoint.push_back(trans(E2, -50, 0));
//				tripoint.push_back(trans(E2, -50, 0)); tripoint.push_back(trans(E3, -50, 0));
//
//				tripoint.push_back(trans(E3, -50, 0)); tripoint.push_back(trans(B3, -50, 0));
//				tripoint.push_back(trans(B3, -50, 0)); tripoint.push_back(trans(A3, -50, 0));
//				tripoint.push_back(trans(A3, -50, 0)); tripoint.push_back(trans(D3, -50, 0));
//				tripoint.push_back(trans(A3, -50, 0)); tripoint.push_back(trans(C3, -50, 0));
//				tripoint.push_back(trans(C3, -50, 0)); tripoint.push_back(trans(F3, -50, 0));
//				tripoint.push_back(trans(C3, -50, 0)); tripoint.push_back(trans(B3, -50, 0));
//
//				tripoint.push_back(trans(B, -50, 0)); tripoint.push_back(trans(A, -50, 0));
//				tripoint.push_back(trans(A, -50, 0)); tripoint.push_back(trans(C, -50, 0));
//				tripoint.push_back(trans(C, -50, 0)); tripoint.push_back(trans(B, -50, 0));
//
//			}
//			n = 0;
//		}
//		break;
//	}
//	case 2:
//	{
//		for (m; m < (120 * 2) / (sqrt(3)*a); m++)   //阶段二：w*i在 sqrt(3)/3 *a——2 * sqrt(3)/3 *a
//		{
//			for (n; n < 170 / a; n++)
//			{
//				point A3 = { m*0.5*x + n * x + 0.5*a ,                           m*y + width * (i - 1) - sqrt(3)*0.333*a + sqrt(3)*0.167*a + 0.5*offset };
//
//				point B3 = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) - 0.25*sqrt(3)*offset ,     m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) - 0.25*offset };
//
//				point C3 = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) + 0.25*sqrt(3)*offset ,       m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) - 0.25*offset };
//
//				point D3 = { m*0.5*x + n * x + 0.5*a ,                           m*y + 0.5*sqrt(3)*a + 0.5*offset };
//
//				point E3 = { m * 0.5 * x + n * x + 0 - 0.25*sqrt(3)*offset,                           m*y + 0 - 0.25*offset };
//
//				point F3 = { m*0.5*x + n * x + a + 0.25*sqrt(3)*offset,                               m*y + 0 - 0.25*offset };
//
//				point A = { m*0.5*x + n * x + 0.5*a ,                           m*y + width * (i - 1) - sqrt(3)*0.333*a + sqrt(3)*0.167*a };
//
//				point A1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*offset,     m*y + width * (i - 1) - sqrt(3)*0.333*a + sqrt(3)*0.167*a + 0.5*offset };
//
//				point A2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*offset,     m*y + width * (i - 1) - sqrt(3)*0.333*a + sqrt(3)*0.167*a + 0.5*offset };
//
//				point B = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) ,     m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };
//
//				point B1 = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) - 0.33*sqrt(3)*offset ,     m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };
//
//				point B2 = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) - 0.167*sqrt(3)*offset ,     m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) - 0.5*offset };
//
//				point C = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) ,       m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };
//
//				point C1 = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) + 0.167*sqrt(3)*offset ,       m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) - 0.5*offset };
//
//				point C2 = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) + 0.33*sqrt(3)*offset,       m*y + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };
//
//				point D1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*offset,                         m*y + 0.5*sqrt(3)*a + 0.5*offset };
//
//				point D2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*offset,                         m*y + 0.5*sqrt(3)*a + 0.5*offset };
//
//				point E1 = { m * 0.5 * x + n * x + 0 - 0.33*sqrt(3)*offset ,                         m*y + 0 };
//
//				point E2 = { m * 0.5 * x + n * x + 0 - 0.167*sqrt(3)*offset ,                         m*y + 0 - 0.5*offset };
//
//				point F1 = { m*0.5*x + n * x + a + 0.33*sqrt(3)*offset,                             m*y + 0 };
//
//				point F2 = { m*0.5*x + n * x + a + 0.167*sqrt(3)*offset,                             m*y + 0 - 0.5*offset };
//
//				tripoint.push_back(trans(E3, -50, 0)); tripoint.push_back(trans(E1, -50, 0));
//				tripoint.push_back(trans(E1, -50, 0)); tripoint.push_back(trans(B1, -50, 0));
//				tripoint.push_back(trans(B1, -50, 0)); tripoint.push_back(trans(A1, -50, 0));
//				tripoint.push_back(trans(A1, -50, 0)); tripoint.push_back(trans(D1, -50, 0));
//				tripoint.push_back(trans(D1, -50, 0)); tripoint.push_back(trans(D2, -50, 0));
//				tripoint.push_back(trans(D2, -50, 0)); tripoint.push_back(trans(A2, -50, 0));
//				tripoint.push_back(trans(A2, -50, 0)); tripoint.push_back(trans(C2, -50, 0));
//				tripoint.push_back(trans(C2, -50, 0)); tripoint.push_back(trans(F1, -50, 0));
//				tripoint.push_back(trans(F1, -50, 0)); tripoint.push_back(trans(F2, -50, 0));
//				tripoint.push_back(trans(F2, -50, 0)); tripoint.push_back(trans(C1, -50, 0));
//				tripoint.push_back(trans(C1, -50, 0)); tripoint.push_back(trans(B2, -50, 0));
//				tripoint.push_back(trans(B2, -50, 0)); tripoint.push_back(trans(E2, -50, 0));
//				tripoint.push_back(trans(E2, -50, 0)); tripoint.push_back(trans(E3, -50, 0));
//
//				tripoint.push_back(trans(E3, -50, 0)); tripoint.push_back(trans(B3, -50, 0));
//				tripoint.push_back(trans(B3, -50, 0)); tripoint.push_back(trans(A3, -50, 0));
//				tripoint.push_back(trans(A3, -50, 0)); tripoint.push_back(trans(D3, -50, 0));
//				tripoint.push_back(trans(A3, -50, 0)); tripoint.push_back(trans(C3, -50, 0));
//				tripoint.push_back(trans(C3, -50, 0)); tripoint.push_back(trans(F3, -50, 0));
//				tripoint.push_back(trans(C3, -50, 0)); tripoint.push_back(trans(B3, -50, 0));
//
//				tripoint.push_back(trans(B, -50, 0)); tripoint.push_back(trans(A, -50, 0));
//				tripoint.push_back(trans(A, -50, 0)); tripoint.push_back(trans(C, -50, 0));
//				tripoint.push_back(trans(C, -50, 0)); tripoint.push_back(trans(B, -50, 0));
//
//			}
//			n = 0;
//		}
//		break;
//	}
//	case 3:
//	{
//		i = i - 2 * sqrt(3)*0.333*a / width;
//		a = a + 0.33*sqrt(3)*offset;
//		x = a + 2 * sqrt(3)*0.33*offset + 0.25*offset;  //x坐标右移单位 这里是外三角形边长
//		y = 0.5*sqrt(3)*a + offset + 0.25*offset;  //y坐标上移单位
//
//		for (m; m < (120 * 2) / (sqrt(3)*a); m++)  //阶段三：和阶段一的区别只是Y坐标全变成负
//		{
//			for (n; n < 170 / a; n++)
//			{
//				point A3 = { m*0.5*x + n * x + 0.5*a ,                         m*y - sqrt(3) * 0.5 * a + width * (i - 1) - 0.5*offset };
//
//				point B3 = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) - 0.25*sqrt(3)*offset,       m*y - 0.5*width * (i - 1) + 0.25*offset };
//
//				point C3 = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) + 0.25*sqrt(3)*offset,   m*y - 0.5*width * (i - 1) + 0.25*offset };
//
//				point D3 = { m*0.5*x + n * x + 0.5*a ,                          m*y - 0.5*sqrt(3)*a - 0.5*offset };
//
//				point E3 = { m * 0.5 * x + n * x + 0 - 0.25*sqrt(3)*offset,                          m*y + 0 + 0.25*offset };
//
//				point F3 = { m*0.5*x + n * x + a + 0.25*sqrt(3)*offset ,                              m*y + 0 + 0.25*offset };
//
//				point A = { m*0.5*x + n * x + 0.5*a ,                         m * y - sqrt(3) * 0.5 * a + width * (i - 1) };
//
//				point A1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*offset ,    m*y - sqrt(3) * 0.5 * a + width * (i - 1) - 0.5*offset };
//
//				point A2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*offset,      m*y - sqrt(3) * 0.5 * a + width * (i - 1) - 0.5*offset };
//
//				point B = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) ,       m*y - 0.5*width * (i - 1) };
//
//				point B1 = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) - 0.33*sqrt(3)*offset ,       m*y - 0.5*width * (i - 1) };
//
//				point B2 = { m*0.5*x + n * x + 0.5*sqrt(3)*width*(i - 1) - 0.167*sqrt(3)*offset ,       m*y - 0.5*width * (i - 1) + 0.5*offset };
//
//				point C = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) ,   m*y - 0.5*width * (i - 1) };
//
//				point C1 = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) + 0.167*sqrt(3)*offset ,   m*y - 0.5*width * (i - 1) + 0.5*offset };
//
//				point C2 = { m*0.5*x + n * x + a - 0.5*sqrt(3)*width*(i - 1) + 0.33*sqrt(3)*offset ,   m*y - 0.5*width * (i - 1) };
//
//				point D1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*offset,                         m*y - 0.5*sqrt(3)*a - 0.5*offset };
//
//				point D2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*offset,                         m*y - 0.5*sqrt(3)*a - 0.5*offset };
//
//				point E1 = { m * 0.5 * x + n * x + 0 - 0.33*sqrt(3)*offset ,                         m*y + 0 };
//
//				point E2 = { m * 0.5 * x + n * x + 0 - 0.167*sqrt(3)*offset ,                         m*y + 0 + 0.5*offset };
//
//				point F1 = { m*0.5*x + n * x + a + 0.33*sqrt(3)*offset,                             m*y + 0 };
//
//				point F2 = { m*0.5*x + n * x + a + 0.167*sqrt(3)*offset,                             m*y + 0 + 0.5*offset };
//
//				tripoint.push_back(trans(E3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(E1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(E1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(B1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(A1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(D1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(D1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(D2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(D2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(A2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(C2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(F1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(F1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(F2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(F2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(C1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(B2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(E2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(E2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(E3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//
//				tripoint.push_back(trans(E3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(B3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(A3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(D3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(A3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(C3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(F3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(C3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//
//				tripoint.push_back(trans(B, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(A, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(C, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//
//				//tripoint.push_back(trans(E3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(E1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(E1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(B1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(A1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(D1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(D1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(D2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(D2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(A2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(C2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(F1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(F1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(F2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(F2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(C1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(B2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(E2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(E2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(E3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//
//				//tripoint.push_back(trans(E3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(B3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(A3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(D3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(A3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(C3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(F3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(C3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//
//				//tripoint.push_back(trans(B, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(A, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				//tripoint.push_back(trans(C, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//
//			}
//			n = 0;
//		}
//		break;
//	}
//	case 4:
//	{
//		i = i - 2 * sqrt(3)*0.333*a / width;
//		a = a + 0.33*sqrt(3)*offset;
//		x = a + 2 * sqrt(3)*0.33*offset + 0.25*offset;  //x坐标右移单位 这里是外三角形边长
//		y = 0.5*sqrt(3)*a + offset + 0.25*offset;  //y坐标上移单位
//
//		for (m; m < (120 * 2) / (sqrt(3)*a); m++)  //阶段四：和阶段二一样
//		{
//			for (n; n < 170 / a; n++)
//			{
//
//				point A3 = { m*0.5*x + n * x + 0.5*a ,                           m*y - width * (i - 1) + sqrt(3)*0.333*a - sqrt(3)*0.167*a - 0.5*offset };
//
//				point B3 = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) - 0.25*sqrt(3)*offset ,     m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) + 0.25*offset };
//
//				point C3 = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) + 0.25*sqrt(3)*offset ,       m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) + 0.25*offset };
//
//				point D3 = { m*0.5*x + n * x + 0.5*a ,                           m*y - 0.5*sqrt(3)*a - 0.5*offset };
//
//				point E3 = { m * 0.5 * x + n * x + 0 - 0.25*sqrt(3)*offset,                           m*y + 0 + 0.25*offset };
//
//				point F3 = { m*0.5*x + n * x + a + 0.25*sqrt(3)*offset,                               m*y + 0 + 0.25*offset };
//
//				point A = { m*0.5*x + n * x + 0.5*a ,                           m*y - width * (i - 1) + sqrt(3)*0.333*a - sqrt(3)*0.167*a };
//
//				point A1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*offset,     m*y - width * (i - 1) + sqrt(3)*0.333*a - sqrt(3)*0.167*a - 0.5*offset };
//
//				point A2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*offset,     m*y - width * (i - 1) + sqrt(3)*0.333*a - sqrt(3)*0.167*a - 0.5*offset };
//
//				point B = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) ,     m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };
//
//				point B1 = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) - 0.33*sqrt(3)*offset ,     m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };
//
//				point B2 = { m*0.5*x + n * x + a - sqrt(3)*0.5*width*(i - 1) - 0.167*sqrt(3)*offset ,     m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) + 0.5*offset };
//
//				point C = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) ,       m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };
//
//				point C1 = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) + 0.167*sqrt(3)*offset ,       m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) + 0.5*offset };
//
//				point C2 = { m*0.5*x + n * x + sqrt(3)*0.5*width*(i - 1) + 0.33*sqrt(3)*offset,       m*y - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };
//
//				point D1 = { m*0.5*x + n * x + 0.5*a - 0.167*sqrt(3)*offset,                         m*y - 0.5*sqrt(3)*a - 0.5*offset };
//
//				point D2 = { m*0.5*x + n * x + 0.5*a + 0.167*sqrt(3)*offset,                         m*y - 0.5*sqrt(3)*a - 0.5*offset };
//
//				point E1 = { m * 0.5 * x + n * x + 0 - 0.33*sqrt(3)*offset ,                         m*y + 0 };
//
//				point E2 = { m * 0.5 * x + n * x + 0 - 0.167*sqrt(3)*offset ,                         m*y + 0 + 0.5*offset };
//
//				point F1 = { m*0.5*x + n * x + a + 0.33*sqrt(3)*offset,                             m*y + 0 };
//
//				point F2 = { m*0.5*x + n * x + a + 0.167*sqrt(3)*offset,                             m*y + 0 + 0.5*offset };
//
//				tripoint.push_back(trans(E3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(E1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(E1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(B1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(A1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(D1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(D1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(D2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(D2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(A2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(C2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(F1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(F1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(F2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(F2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C1, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(C1, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(B2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(E2, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(E2, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(E3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//
//				tripoint.push_back(trans(E3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(B3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(A3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(D3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(A3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(C3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(F3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(C3, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B3, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//
//				tripoint.push_back(trans(B, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(A, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(A, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(C, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//				tripoint.push_back(trans(C, -50 - 0.167*sqrt(3)*offset, -0.5*offset)); tripoint.push_back(trans(B, -50 - 0.167*sqrt(3)*offset, -0.5*offset));
//
//			}
//			n = 0;
//		}
//		break;
//	}
//
//	default: printf("有问题");
//		break;
//	}
//
//}

void intersection()           //在这里面调节要不要变化
{
	for (int i = 0; i < model.size(); i++)
	{
		tripoint.clear(); //把上一层的变化线清理掉
		interpoint.clear(); //把上一层的交点清理掉
		polypoint.clear();//把上一层的轮廓线清理掉
		triangle(i); //产生第i层的变化线存在tripoint        0：不变三角，48不变六边形， i变化线

		polypoint = model[i];//把第i层的轮廓线给polypoint

		for (int m = 0; m < tripoint.size(); m += 2)  //网格的每一条线段，首尾相连线段端点重复，所以跳两个
		{
			interpoint.clear(); //上一条线段的交点清理掉
			for (int n = 0; n < polypoint.size(); n++) //多边形的每一条线段，首尾相连线段端点不重复，所以跳一个
			{
				inter(tripoint[m], tripoint[m + 1], polypoint[n], polypoint[(n + 1 + polypoint.size()) % polypoint.size()]);
			}
			int a = InOrOutPolygon(tripoint[m]);
			int b = InOrOutPolygon(tripoint[m + 1]);  //判断是不是内点,0:外 1:内

			if (interpoint.size() == 0)   // 无交点
			{
				if ((a == 0) && (b == 0)) //无内点，舍
				{
					continue;
				}
				if ((a == 1) || (b == 1)) //存在内点，取
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(tripoint[m + 1]);
					continue;
				}
			}

			if (interpoint.size() == 1)   //有一交点
			{
				if ((a == 0) && (b == 0))   //有一交点且无内点，舍 （相切）
				{
					continue;
				}
				if ((a == 1) && (b == 0))//有一交点且有一内点，取内点到交点 （相交）
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(interpoint[0]);
					continue;
				}
				if ((a == 0) && (b == 1))//有一交点且有一内点，取内点到交点 （相交）
				{
					new_tripoint.push_back(interpoint[0]);
					new_tripoint.push_back(tripoint[m + 1]);
					continue;
				}
				if ((a == 1) && (b == 1)) //有一交点且有两内点，取 （内部包含）
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(tripoint[m + 1]);
					continue;
				}

			}

			if (interpoint.size() == 2)   //有两交点
			{
				if ((a == 0) && (b == 0))   //有两交点且无内点，取两交点
				{
					new_tripoint.push_back(interpoint[0]);
					new_tripoint.push_back(interpoint[1]);
					continue;
				}
				if ((a == 1) && (b == 0))//有两交点且一个内点 两两相连，注意顺序 
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(mindistance(tripoint[m], interpoint[0], interpoint[1]));
					continue;
				}
				if ((a == 0) && (b == 1))//有两交点且一个内点 两两相连，注意顺序
				{
					new_tripoint.push_back(tripoint[m + 1]);
					new_tripoint.push_back(mindistance(tripoint[m + 1], interpoint[0], interpoint[1]));
					continue;
				}
				if ((a == 1) && (b == 1))//有两交点且两个内点 两两相连，注意顺序 
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(mindistance(tripoint[m], interpoint[0], interpoint[1]));
					new_tripoint.push_back(mindistance(tripoint[m + 1], interpoint[0], interpoint[1]));
					new_tripoint.push_back(tripoint[m + 1]);
					continue;
				}

			}

		}

		modelfill.push_back(new_tripoint);
		new_tripoint.clear();
		printf("%.2lf%%\r", i * 100.0 / model.size());
	}


}

void writegcode()
{
	FILE* fp;

	errno_t err;     //判断此文件流是否存在 存在返回1

	err = fopen_s(&fp, "test gcode.txt", "a"); //若return 1 , 则将指向这个文件的文件流给fp1

	fprintf(fp, "M104 S200\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M109 S200\n");
	fprintf(fp, "M82;set extruder to absolute mode\n");
	fprintf(fp, "G28;move X/Y to min endstops\n");
	fprintf(fp, "G1 Z15.0 F6000 ;move the platform down 15mm\n");
	fprintf(fp, "G92 E0                  ;zero the extruded length\n");
	fprintf(fp, "G1 F200 E3             ;extrude 3mm of feed stock\n");
	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G1 F1500 E-6.5\n");

	fprintf(fp, ";LAYER_COUNT: %d\n", modelfill.size() + 1);
	double E = 0;

	fprintf(fp, ";LAYER:%d\n", 0);
	//for (int j = 0; j < 25; j++) //大底座
	//{
	//	fprintf(fp, "G0 F300 X%.3f Y%.3f Z%.3f\n", j * 4 + 10.00, 10.00, 0.500);
	//	fprintf(fp, "G1 F1000 X%.3f Y%.3f E%.3f\n", j * 4 + 10.00, 100.00, E += distance({ j * 4 + 10.00, 10.00 }, { j * 4 + 10.00, 100.00 })*0.3);
	//	fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 12.00, 100.00, E += distance({ j * 4 + 10.00, 100.00 }, { j * 4 + 12.00, 100.00 })*0.3);
	//	fprintf(fp, "G1 F1000 X%.3f Y%.3f E%.3f\n", j * 4 + 12.00, 10.00, E += distance({ j * 4 + 12.00, 100.00 }, { j * 4 + 12.00, 10.00 })*0.3);
	//	fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 14.00, 10.00, E += distance({ j * 4 + 12.00, 10.00 }, { j * 4 + 14.00, 10.00 })*0.3);
	//}

	for (int j = 0; j < 15; j++) //小底座
	{
		fprintf(fp, "G0 F300 X%.3f Y%.3f Z%.3f\n", j * 4 + 30.00, 30.00, 0.500);
		fprintf(fp, "G1 F1000 X%.3f Y%.3f E%.3f\n", j * 4 + 30.00, 100.00, E += distance({ j * 4 + 30.00, 30.00 }, { j * 4 + 30.00, 100.00 })*0.3);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 32.00, 100.00, E += distance({ j * 4 + 30.00, 100.00 }, { j * 4 + 32.00, 100.00 })*0.3);
		fprintf(fp, "G1 F1000 X%.3f Y%.3f E%.3f\n", j * 4 + 32.00, 30.00, E += distance({ j * 4 + 32.00, 100.00 }, { j * 4 + 32.00, 30.00 })*0.3);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 34.00, 30.00, E += distance({ j * 4 + 32.00, 30.00 }, { j * 4 + 34.00, 30.00 })*0.3);
	}

	for (int i = 1; i < model.size(); i++) //每一层 i代表层数      z层厚 t挤出率
	{
		fprintf(fp, ";LAYER:%d\n", i);

		//fprintf(fp, "G0 F2000 Z%.3f\n", 0.500 + i * z);//去掉轮廓线时用到，确定Z

		fprintf(fp, "G0 F2000 X%.3f Y%.3f Z%.3f\n", model[i][0].x, model[i][0].y, 0.500 + i * z);
		fprintf(fp, ";TYPE:OUTLINE\n");
		fprintf(fp, "G1 F1200 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][0], model[i][1])*t);
		for (int j = 2; j < model[i].size(); j++)
		{
			fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][j].x, model[i][j].y, E += distance(model[i][j - 1], model[i][j])*t);
		}
		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][0].x, model[i][0].y, E += distance(model[i][model[i].size() - 1], model[i][0])*t); //画一个圈，要回原点
		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][model[i].size() - 1], model[i][0])*0.01);//回到原点再多出去一点，有利于不翘边

		fprintf(fp, ";TYPE:FILL\n");
		for (int k = 1; k < modelfill[i].size(); k += 2)
		{
			fprintf(fp, "G0 F2000 X%.3f Y%.3f\n", modelfill[i][k - 1].x, modelfill[i][k - 1].y);
			fprintf(fp, "G1 F1200 X%.3f Y%.3f E%.5f\n", modelfill[i][k].x, modelfill[i][k].y, E += distance(modelfill[i][k - 1], modelfill[i][k])*t);
		}
	}

	fprintf(fp, "M107\n");
	fprintf(fp, "M104 S0                     ;extruder heater off\n");
	fprintf(fp, "G92 E1                      ;relative positioning\n");
	fprintf(fp, "G1 E-1 F300                 ;retract the filament a bit before lifting the nozzle, to release some of the pressure\n");
	fprintf(fp, "G1 X120 Y120 Z120 F9000     ;move to max so the object can take out\n");
	fprintf(fp, "M84                         ;steppers off\n");
	fprintf(fp, "M82                         ;absolute extrusion mode\n");
	fprintf(fp, "M104 S0\n");
	fprintf(fp, ";End GCode\n");

	fclose(fp);
}

void writpointcloud()
{
	FILE* fp;

	errno_t err;     //判断此文件流是否存在 存在返回1

	err = fopen_s(&fp, "test pointcloud.txt", "a"); //若return 1 , 则将指向这个文件的文件流给fp1

	double E = 0;

	for (int i = 0; i < model.size(); i++) //每一层 i代表层数      z层厚 t挤出率
	{

		//for (int j = 1; j < model[i].size(); j += 2)
		//{
		//	fprintf(fp, "%.3f %.3f %.3f\t", model[i][j-1].x, model[i][j-1].y, i*0.2);
		//	fprintf(fp, "%.3f %.3f %.3f\n", model[i][j].x, model[i][j].y, i*0.2);
		//}

		for (int k = 1; k < modelfill[i].size(); k += 2)
		{
			fprintf(fp, "%.3f %.3f %.3f\t", modelfill[i][k - 1].x, modelfill[i][k - 1].y, i*0.2);
			fprintf(fp, "%.3f %.3f %.3f\n", modelfill[i][k].x, modelfill[i][k].y,  i*0.2);
		}
	}

	fclose(fp);
}

void InitGL(GLvoid)
{
	glShadeModel(GL_SMOOTH);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_COLOR_MATERIAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
}

void lines() //画线的函数


{
	//glColor3f(1.0, 0.0, 0.0);
	//glBegin(GL_LINES);
	//for (int j = 0; j < 50; j++)
	//{
	//	for (int i = 0; i < tripoint.size(); i += 2)
	//	{
	//		glVertex3f(tripoint[i].x, tripoint[i].y, 0);  glVertex3f(tripoint[i + 1].x, tripoint[i + 1].y, 0);
	//	}
	//}
	//glEnd();
	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);  glVertex3f(0, 0, 0); glVertex3f(100, 0, 0);
	glColor3f(0.0, 1.0, 0.0); glVertex3f(0, 0, 0); glVertex3f(0, 100, 0);
	glColor3f(0.0, 0.0, 1.0); glVertex3f(0, 0, 0); glVertex3f(0, 0, 100);
	glEnd();
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	for (int i = 1; i < modelfill.size(); i++)
	{
		for (int j = 1; j < modelfill[i].size(); j += 2)
		{
			glVertex3f(modelfill[i][j-1].x, modelfill[i][j-1].y,i*0.2); glVertex3f(modelfill[i][j].x, modelfill[i][j].y, i*0.2);
		}
	}
	glEnd();

}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glPushMatrix();

	gluLookAt(60 + xs, 60 + ys, 200 + zs, 60 + sx, 60 + sy, 60 + sz, 0, 1, 0);
	//glTranslatef(0.0f, 0.0f, -1.0f);	//平移
	//glRotatef(rquad, 1.0f, 0.0f, 0.0f);	//旋转一个角度
	lines();

	glColor3f(1.0, 0.0, 0.0);
	glPointSize(5);

	glBegin(GL_POINTS);
	//glVertex3f(0, 0, 0);
	//glVertex3f(-50, 0, 1);
	glVertex3f(60,60,60);
	glEnd();
	glPopMatrix();

	//rquad -= 0.15f;	//修改立方体的旋转角度
	glutSwapBuffers();
}

void reshape(int width, int height)
{
	if (height == 0)
		height = 1;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(90.0f, (GLfloat)width / (GLfloat)height, 0.0f, 100.0f);
	gluPerspective(90.0f, (GLfloat)width / (GLfloat)height, 0.0f, 100.0f);
	//if (width <= height)
	//	glOrtho(-2.0, 2.0, -2.0*(GLfloat)height / (GLfloat)width, 2.0*(GLfloat)height / (GLfloat)width, 1.0, 20.0);
	//else
	//	glOrtho(-2.0*(GLfloat)width / (GLfloat)height, 2.0*(GLfloat)width / (GLfloat)height, -2.0, 2.0, 1.0, 20.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void myMouse(int button, int state, int x, int y)
{
	//if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
	//	mousetate = 1;
	//	Oldx = x;
	//	Oldy = y;
	//}
	//if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
	//	mousetate = 0;
	//滚轮事件
	if (state == GLUT_UP && button == 3) {

		zs -= 1;
	}
	if (state == GLUT_UP && button == 4) {

		zs += 1;
	}
	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'w':
		sy += 10;
		glutPostRedisplay();
		break;
	case 's':
		sy -= 10;
		glutPostRedisplay();
		break;
	case 'a':
		sx -= 10;
		glutPostRedisplay();
		break;
	case 'd':
		sx += 10;
		glutPostRedisplay();
		break;
	case 'q':
		sz += 10;
		glutPostRedisplay();
		break;
	case 'e':
		sz -= 10;
		glutPostRedisplay();
		break;
	default:
		break;
	}
}

void SpecialKey(GLint key, GLint x, GLint y)
{
	if (key == GLUT_KEY_UP)
	{
		ys += 10;
	}
	if (key == GLUT_KEY_LEFT)
	{
		xs -= 10;
	}
	if (key == GLUT_KEY_DOWN)
	{
		ys -= 10;
	}
	if (key == GLUT_KEY_RIGHT)
	{
		xs += 10;
	}
	glutPostRedisplay();
}

void main(int argc, char** argv)
{
	//************************************************************************切片
	printf("slicing...\n");
	readfile(file_2);
	BoundingBox();
	findIntersect();  //前三个函数切片，产生model
	printf("slice complete,layer count: %d\n", model.size());

	//************************************************************************取交
	printf("path planning...\n");
	intersection();
	printf("path planning complete\n"); //modelfill完成

	//************************************************************************编写Gcode
	printf("Gcode writing...\n");
	//writegcode();
	writpointcloud();
	printf("Gcode writing complete,file save as \"test gcode.txt\"\n");

	//************************************************************************预览
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(600, 600);
	glutCreateWindow("Hello Cube");
	InitGL();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutIdleFunc(display);
	glutMouseFunc(myMouse);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(&SpecialKey);
	glutMainLoop();

}


