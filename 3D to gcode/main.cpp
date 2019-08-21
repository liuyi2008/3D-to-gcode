//#include<stdio.h>
//#include<iostream>
//#include <vector>
//#include<math.h>
//#include<windows.h>
//#include<GL/glut.h>

#include"slice.h"

typedef unsigned char boolean;

#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)

//struct point { double x, y; }p1;

point p1;

vector<point>polypoint; //一层的轮廓线坐标
vector<point>tripoint;  //一层的变化线 线段坐标 成对
vector<point>new_tripoint;  //新的一层的变化线 线段坐标 成对
vector<point>interpoint; //交点坐标

extern const string file_1;  //slice
extern const string file_2; 

vector<vector<point> >model; //slice  整个模型的，分层的轮廓线坐标，一维层数，二维此层的罗廓线点坐标
vector<vector<point> >modelfill(model); //完成的model填充线，一维层数，二维此层的变化线 线段坐标

//**************************************************************判断线段与线段的关系函数********************************************************
//计算交叉乘积(P1-P0)x(P2-P0)
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

//****************************************************************判断线段与线段的关系函数************************************************************

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
	return sqrt( pow((a.x - b.x), 2) + pow((a.y - b.y), 2) );
}

void triangle(int i)  //第i层，从1开始，由i和w一起决定点A、B、C的位置
{
	float Thickness = 0.2; // 一层的厚度

	float a = 100.0f;//单位距离 这里来控制缩放

	float width = 0.2; //喷头的宽度

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
				//tripoint.push_back(trans(B)); tripoint.push_back(trans(E));

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
				//tripoint.push_back(trans(B)); tripoint.push_back(trans(E));

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
				//tripoint.push_back(trans(B)); tripoint.push_back(trans(E));

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
				//tripoint.push_back(trans(B)); tripoint.push_back(trans(E));

			}
			n = 0;
		}
		break;
	}

	default: printf("有问题");
		break;
	}

}

//****************************************************************画出来测试用************************************************************

void init(int argc, char** argv)
{
	glutInit(&argc, argv);  //I初始化 GLUT.
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);  //设置显示模式：单个缓存和使用RGB模型
	glutInitWindowPosition(50, 100);  //设置窗口的顶部和左边位置
	glutInitWindowSize(400, 300);  //设置窗口的高度和宽度
	glutCreateWindow("Scan Program");

	glClearColor(1.0, 1.0, 1.0, 0); //窗口背景颜色为白色
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, 600, 0, 450);
}

void myDisplay(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glPointSize(1);
	glBegin(GL_LINES);

	for (int i = 0; i < new_tripoint.size(); i++)
	{

		glColor3f(1.0f, 0.0f, 0.0f); //线是红色的
		glVertex2f(new_tripoint[i].x, new_tripoint[i].y);
		//printf("i=%d,x=%f,y=%f\n",i , tripoint[i].x, tripoint[i].y);		

	}

	glEnd();
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < polypoint.size(); i++)
	{

		glColor3f(0.0f, 0.0f, 1.0f); //线是红色的
		glVertex2f(polypoint[i].x, polypoint[i].y);
		//printf("i=%d,x=%f,y=%f\n",i , tripoint[i].x, tripoint[i].y);		

	}
	glEnd();
	glFlush();
}

//****************************************************************画出来测试用************************************************************

void main(int argc, char** argv)
{
	printf("slicing...\n");
	readfile(file_2);
	BoundingBox();
	findIntersect();  //前三个函数切片，产生model

	//for (int i = 0; i < model.size(); i++)
	//{
	//	for (int j = 0; j < model[i].size(); j++)
	//	{
	//		printf("x = %d, y = %f,z = %f\n", i, model[i][j].x, model[i][j].y);
	//	}
	//}

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

			p1 = { tripoint[i].x*0.2,tripoint[i].y*0.2};
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


	/*init(argc, argv);
	glutDisplayFunc(myDisplay);
	glutMainLoop();*/

	//*******************************************以下是编写Gcode*******************************************

	printf("Gcode writing...\n");

	FILE* fp;

	errno_t err;     //判断此文件流是否存在 存在返回1

	err = fopen_s(&fp, "test gcode.txt", "a"); //若return 1 , 则将指向这个文件的文件流给fp1

	fprintf(fp, "M104 S190\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M109 S190\n");
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
	for (int j = 0; j < 18; j++) //底座
	{
		fprintf(fp, "G0 F300 X%.3f Y%.3f Z%.3f\n", j * 4 + 30.00, 30.00,  0.500);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 30.00, 100.00, E += distance({ j * 4 + 30.00, 30.00 }, { j * 4 + 30.00, 100.00 })*0.3);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 32.00, 100.00, E += distance({ j * 4 + 30.00, 100.00 }, { j * 4 + 32.00, 100.00 })*0.3);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 32.00, 30.00,  E += distance({ j * 4 + 32.00, 100.00 }, { j * 4 + 32.00, 30.00 })*0.3);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 34.00, 30.00,  E += distance({ j * 4 + 32.00, 30.00 }, { j * 4 + 34.00, 30.00 })*0.3);

	}

	//for (int i = 0; i < 1; i++) //第一层打的粗一些，相当于底座
	//{
	//	fprintf(fp, ";LAYER:%d\n", i);

	//	fprintf(fp, "G0 F1000 X%.3f Y%.3f Z%.3f\n", 0.00, 0.00, 0.400 + i * 0.200);  //原点
	//	fprintf(fp, "G1 F1000 X%.3f Y%.3f E%.3f\n", model[i][0].x, model[i][0].y, E += distance({0.0}, model[i][0])*0.6);
	//	fprintf(fp, ";TYPE:OUTLINE\n");
	//	for (int j = 1; j < model[i].size(); j++)
	//	{
	//		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][j].x, model[i][j].y, E += distance(model[i][j - 1], model[i][j])*0.6);
	//	}
	//	fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][0].x, model[i][0].y, E += distance(model[i][model[i].size() - 1], model[i][0])*0.6); //画一个圈，要回原点
	//	fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][model[i].size() - 1], model[i][0])*0.01);

	//	fprintf(fp, ";TYPE:FILL\n");
	//	for (int k = 1; k < modelfill[i].size(); k += 2)
	//	{
	//		fprintf(fp, "G0 F2600 X%.3f Y%.3f\n", modelfill[i][k - 1].x, modelfill[i][k - 1].y);
	//		fprintf(fp, "G1 F1000 X%.3f Y%.3f E%.5f\n", modelfill[i][k].x, modelfill[i][k].y, E += distance(modelfill[i][k - 1], modelfill[i][k])*0.1);
	//	}
	//}

	for (int i = 1; i < model.size(); i++) //每一层 i代表层数
	{
		fprintf(fp, ";LAYER:%d\n", i);	
		fprintf(fp, "G0 F2000 X%.3f Y%.3f Z%.3f\n", model[i][0].x, model[i][0].y, 0.500 + i * 0.300);
		fprintf(fp, ";TYPE:OUTLINE\n");
		fprintf(fp, "G1 F1000 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][0], model[i][1])*0.1);
		for (int j = 2; j < model[i].size(); j++)
		{
			fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n",model[i][j].x, model[i][j].y, E += distance(model[i][j-1], model[i][j])*0.1);
		}
		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][0].x, model[i][0].y, E += distance(model[i][model[i].size()-1], model[i][0])*0.1); //画一个圈，要回原点

		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][model[i].size()-1], model[i][0])*0.01);
		fprintf(fp, ";TYPE:FILL\n");
		for (int k = 1; k < modelfill[i].size(); k+=2)
		{
			fprintf(fp, "G0 F2000 X%.3f Y%.3f\n", modelfill[i][k-1].x, modelfill[i][k-1].y);
			fprintf(fp, "G1 F1000 X%.3f Y%.3f E%.5f\n", modelfill[i][k].x, modelfill[i][k].y, E += distance(modelfill[i][k-1], modelfill[i][k])*0.1);
		}
	}

	fprintf(fp, "M107\n");
	fprintf(fp, "M104 S0                     ;extruder heater off\n");
	fprintf(fp, "M140 S0                     ;heated bed heater off (if you have it)\n");
	fprintf(fp, "G92 E1                                    ;relative positioning\n");
	fprintf(fp, "G1 E-1 F300                            ;retract the filament a bit before lifting the nozzle, to release some of the pressure\n");
	fprintf(fp, "G1 Z+0.5 E-5 X-20 Y-20 F9000 ;move Z up a bit and retract filament even more\n");
	fprintf(fp, "G28 X0 Y0                              ;move X/Y to min endstops, so the head is out of the way\n");
	fprintf(fp, "M84                         ;steppers off\n");
	fprintf(fp, "M82 ;absolute extrusion mode\n");
	fprintf(fp, "M104 S0\n");
	fprintf(fp, ";End GCode\n");

	fclose(fp);

	printf("Gcode writing complete,file save as \"test gcode.txt\"\n");
}
