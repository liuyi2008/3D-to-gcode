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

}p,p1;

typedef unsigned char boolean;

#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)

vector<point>coord;
vector<point>polypoint; //一层的轮廓线坐标
vector<point>tripoint;  //一层的变化线 线段坐标 成对
vector<point>new_tripoint;  //新的一层的变化线 线段坐标 成对
vector<point>interpoint; //交点坐标
vector<vector<point> >model; //slice  整个模型的，分层的轮廓线坐标，一维层数，二维此层的罗廓线点坐标
vector<vector<point> >modelfill(model); //完成的model填充线，一维层数，二维此层的变化线 线段坐标

const string file_1 = "test.obj";  //slice
const string file_2 = "buddha_head.obj";
const string file_3 = "wawa.obj";
const string file_4 = "opener.obj";
const string file_5 = "test2.obj";
const string file_6 = "bunnyr.obj";
int currentfile = 1;

const float z = 0.30;   //层厚

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
	for (double i = Bmin_y; i < Bmax_y; ) {
		MyMesh::Point pt;
		MyMesh::Point pilist[1000];
		MyMesh::Normal vf(0, 1, 0);
		int pnum = 0;
		float Xport = 0;
		pt.data()[0] = 0; pt.data()[1] = i; pt.data()[2] = 0;
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
					p = { pilist[j].data()[0] ,pilist[j].data()[2] };
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

void inter(point u1, point u2, point v1, point v2)  //主函数，有交点存交点，没交点无操作
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

point trans(point p)  //坐标平移  这里似乎不应该缩放，只需要平移
{
	int b = -250;  //左右平移
	int c = 0;  //上下平移
	p = { p.x + b , p.y + c };
	return p;
}

point trans2(point p)//调整轮廓线
{
	int a = 1;  //缩放
	int b = -250;  //左右平移
	int c = 0;  //上下平移
	p = { a * p.x + b , a * p.y + c };
	return p;
}

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

void triangle(int i)  //第i层，从1开始，由i和w一起决定点A、B、C的位置
{
	float Thickness = 0.2; // 一层的厚度

	float a = 100.0f;//单位距离 这里来控制缩放

	float width = 0.4; //喷头的宽度,或者变化的速率

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
		for (m; m < 10; m++)  //阶段一：w*i 在 0——sqrt(3)/3 *a
		{
			for (n; n < 10; n++)
			{
				point A = { m*0.5*a + n * a + 0.5*a ,                         m*0.5*sqrt(3)*a + sqrt(3) * 0.5 * a - width * (i - 1) };

				point B = { m*0.5*a + n * a + 0.5*sqrt(3)*width*(i - 1) ,       m*0.5*sqrt(3)*a + 0.5*width * (i - 1) };

				point C = { m*0.5*a + n * a + a - 0.5*sqrt(3)*width*(i - 1) ,   m*0.5*sqrt(3)*a + 0.5*width * (i - 1) };

				point D = { m*0.5*a + n * a + 0.5*a ,                         m*0.5*sqrt(3)*a + 0.5*sqrt(3)*a };

				point E = { m * 0.5 * a + n * a + 0 ,                         m*0.5*sqrt(3)*a + 0 };

				point F = { m*0.5*a + n * a + a ,                             m*0.5*sqrt(3)*a + 0 };

				tripoint.push_back(trans(E)); tripoint.push_back(trans(B));
				tripoint.push_back(trans(B)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(F));
				//tripoint.push_back(trans(F)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(D));
				//tripoint.push_back(trans(D)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(B));
				tripoint.push_back(trans(B)); tripoint.push_back(trans(E));

			}
			n = 0;
		}
		break;
	}
	case 2:
	{
		for (m; m < 10; m++)   //阶段二：w*i在 sqrt(3)/3 *a——2 * sqrt(3)/3 *a
		{
			for (n; n < 10; n++)
			{
				point A = { m*0.5*a + n * a + 0.5*a ,                           m*0.5*sqrt(3)*a + width * (i - 1) - sqrt(3)*0.333*a + sqrt(3)*0.167*a };

				point B = { m*0.5*a + n * a + a - sqrt(3)*0.5*width*(i - 1) ,     m*0.5*sqrt(3)*a + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };

				point C = { m*0.5*a + n * a + sqrt(3)*0.5*width*(i - 1) ,       m*0.5*sqrt(3)*a + sqrt(3) * 0.333 * a - 0.5*width * (i - 1) };

				point D = { m*0.5*a + n * a + 0.5*a ,                           m*0.5*sqrt(3)*a + 0.5*sqrt(3)*a };

				point E = { m * 0.5 * a + n * a + 0 ,                           m*0.5*sqrt(3)*a + 0 };

				point F = { m*0.5*a + n * a + a ,                               m*0.5*sqrt(3)*a + 0 };

				tripoint.push_back(trans(E)); tripoint.push_back(trans(B));
				tripoint.push_back(trans(B)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(F));
				//tripoint.push_back(trans(F)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(D));
				//tripoint.push_back(trans(D)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(B));
				tripoint.push_back(trans(B)); tripoint.push_back(trans(E));

			}
			n = 0;
		}
		break;
	}
	case 3:
	{
		i = i - 2 * sqrt(3)*0.333*a / width;
		for (m; m < 10; m++)  //阶段三：和阶段一的区别只是Y坐标全变成负
		{
			for (n; n < 10; n++)
			{
				point A = { m*0.5*a + n * a + 0.5*a ,                         m*0.5*sqrt(3)*a - sqrt(3) * 0.5 * a + width * (i - 1) };

				point B = { m*0.5*a + n * a + 0.5*sqrt(3)*width*(i - 1) ,       m*0.5*sqrt(3)*a - 0.5*width * (i - 1) };

				point C = { m*0.5*a + n * a + a - 0.5*sqrt(3)*width*(i - 1) ,   m*0.5*sqrt(3)*a - 0.5*width * (i - 1) };

				point D = { m*0.5*a + n * a + 0.5*a ,                          m*0.5*sqrt(3)*a - 0.5*sqrt(3)*a };

				point E = { m * 0.5 * a + n * a + 0 ,                          m*0.5*sqrt(3)*a + 0 };

				point F = { m*0.5*a + n * a + a ,                              m*0.5*sqrt(3)*a + 0 };

				tripoint.push_back(trans(E)); tripoint.push_back(trans(B));
				tripoint.push_back(trans(B)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(F));
				//tripoint.push_back(trans(F)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(D));
				//tripoint.push_back(trans(D)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(B));
				tripoint.push_back(trans(B)); tripoint.push_back(trans(E));

			}
			n = 0;
		}
		break;
	}
	case 4:
	{
		i = i - 2 * sqrt(3)*0.333*a / width;
		for (m; m < 10; m++)  //阶段四：和阶段二一样
		{
			for (n; n < 10; n++)
			{
				point A = { m*0.5*a + n * a + 0.5*a ,                           m*0.5*sqrt(3)*a - width * (i - 1) + sqrt(3)*0.333*a - sqrt(3)*0.167*a };

				point B = { m*0.5*a + n * a + a - sqrt(3)*0.5*width*(i - 1) ,     m*0.5*sqrt(3)*a - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };

				point C = { m*0.5*a + n * a + sqrt(3)*0.5*width*(i - 1) ,       m*0.5*sqrt(3)*a - sqrt(3) * 0.333 * a + 0.5*width * (i - 1) };

				point D = { m*0.5*a + n * a + 0.5*a ,                           m*0.5*sqrt(3)*a - 0.5*sqrt(3)*a };

				point E = { m * 0.5 * a + n * a + 0 ,                           m*0.5*sqrt(3)*a + 0 };

				point F = { m*0.5*a + n * a + a ,                               m*0.5*sqrt(3)*a + 0 };

				tripoint.push_back(trans(E)); tripoint.push_back(trans(B));
				tripoint.push_back(trans(B)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(F));
				//tripoint.push_back(trans(F)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(D));
				//tripoint.push_back(trans(D)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(B));
				tripoint.push_back(trans(B)); tripoint.push_back(trans(E));

			}
			n = 0;
		}
		break;
	}

	default: printf("有问题");
		break;
	}

}


void main(int argc, char** argv)
{
	printf("slicing...\n");
	readfile(file_6);
	BoundingBox();
	findIntersect();  //前三个函数切片，产生model
	printf("slice complete,layer count: %d\n", model.size());


	printf("path planning...\n");
	for (int i = 0; i < model.size(); i++)
	{
		tripoint.clear(); //把上一层的变化线清理掉
		interpoint.clear(); //把上一层的交点清理掉
		polypoint.clear();//把上一层的轮廓线清理掉
		triangle(i); //产生第i层的变化线存在tripoint
		for (int i = 0; i < tripoint.size(); i++) //对变化线缩放
		{

			p1 = { tripoint[i].x*0.2,tripoint[i].y*0.2 };
			tripoint[i] = p1;

		}
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
		printf("%.4lf%%\r", i * 100.0 / model.size());
	}

	printf("path planning complete\n"); //modelfill完成

	//*******************************************以下是编写Gcode*******************************************

	printf("Gcode writing...\n");

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
	//for (int j = 0; j < 18; j++) //大底座
	//{
	//	fprintf(fp, "G0 F300 X%.3f Y%.3f Z%.3f\n", j * 4 + 10.00, 10.00, 0.500);
	//	fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 10.00, 100.00, E += distance({ j * 4 + 10.00, 10.00 }, { j * 4 + 10.00, 100.00 })*0.3);
	//	fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 12.00, 100.00, E += distance({ j * 4 + 10.00, 100.00 }, { j * 4 + 12.00, 100.00 })*0.3);
	//	fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 12.00, 10.00, E += distance({ j * 4 + 12.00, 100.00 }, { j * 4 + 12.00, 10.00 })*0.3);
	//	fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 14.00, 10.00, E += distance({ j * 4 + 12.00, 10.00 }, { j * 4 + 14.00, 10.00 })*0.3);

	//}
	for (int j = 0; j < 18; j++) //小底座
	{
		fprintf(fp, "G0 F300 X%.3f Y%.3f Z%.3f\n", j * 4 + 30.00, 30.00, 0.500);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 30.00, 100.00, E += distance({ j * 4 + 30.00, 30.00 }, { j * 4 + 30.00, 100.00 })*0.3);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 32.00, 100.00, E += distance({ j * 4 + 30.00, 100.00 }, { j * 4 + 32.00, 100.00 })*0.3);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 32.00, 30.00, E += distance({ j * 4 + 32.00, 100.00 }, { j * 4 + 32.00, 30.00 })*0.3);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 34.00, 30.00, E += distance({ j * 4 + 32.00, 30.00 }, { j * 4 + 34.00, 30.00 })*0.3);

	}

	for (int i = 1; i < model.size(); i++) //每一层 i代表层数
	{
		fprintf(fp, ";LAYER:%d\n", i);
		//fprintf(fp, "G0 F2000 X%.3f Y%.3f Z%.3f\n", modelfill[i][0].x, modelfill[i][0].y, 0.500 + i * z);  去掉轮廓线时用到了

		fprintf(fp, "G0 F2000 X%.3f Y%.3f Z%.3f\n", model[i][0].x, model[i][0].y, 0.500 + i * z);
		fprintf(fp, ";TYPE:OUTLINE\n");
		fprintf(fp, "G1 F1200 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][0], model[i][1])*0.1);
		for (int j = 2; j < model[i].size(); j++)
		{
			fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][j].x, model[i][j].y, E += distance(model[i][j - 1], model[i][j])*0.1);
		}
		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][0].x, model[i][0].y, E += distance(model[i][model[i].size() - 1], model[i][0])*0.1); //画一个圈，要回原点
		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][model[i].size() - 1], model[i][0])*0.01);

		fprintf(fp, ";TYPE:FILL\n");
		for (int k = 1; k < modelfill[i].size(); k += 2)
		{
			fprintf(fp, "G0 F2000 X%.3f Y%.3f\n", modelfill[i][k - 1].x, modelfill[i][k - 1].y);
			fprintf(fp, "G1 F1200 X%.3f Y%.3f E%.5f\n", modelfill[i][k].x, modelfill[i][k].y, E += distance(modelfill[i][k - 1], modelfill[i][k])*0.1);
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

	printf("Gcode writing complete,file save as \"test gcode.txt\"\n");
}
