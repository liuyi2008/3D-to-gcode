#include<stdio.h>
#include<iostream>
#include <vector>
#include<math.h>
#include<windows.h>
#include<GL/glut.h>



using std::vector;
using namespace std;



#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)

struct point { double x, y; }p1;

vector<point>polypoint = { {250,50},{550,150},{550,400},{250,250},{100,350},{100,100},{120,30} }; //轮廓线坐标集
vector<point>tripoint;  //变化线 线段坐标 成对
vector<point>new_tripoint;  //新变化线 线段坐标 成对
vector<point>interpoint; //交点坐标

//extern const string file_1;  //slice
//vector<vector<point> >model; //slice



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
		bool cond1 = (polypoint[i].x <= x0) && (x0 < polypoint[(i + 1 + n) % n].x);
		bool cond2 = (polypoint[(i + 1 + n) % n].x <= x0) && (x0 < polypoint[i].x);
		bool above = (y0 < slope * (x0 - polypoint[i].x) + polypoint[i].y);
		if ((cond1 || cond2) && above) crossings++;
	}
	return (crossings % 2 != 0);    //返 回 值:  0:外 1:内
}

point trans(point p)  //坐标平移  这里似乎不应该缩放，只需要平移
{
	int b = -250;  //左右平移
	int c = 0;  //上下平移
	p = {p.x + b , p.y + c };
	return p;
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
		for (m; m < 5; m++)  //阶段一：w*i 在 0——sqrt(3)/3 *a
		{
			for (n; n < 8; n++)
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
				tripoint.push_back(trans(F)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(D));
				tripoint.push_back(trans(D)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(B));
				tripoint.push_back(trans(B)); tripoint.push_back(trans(E));

			}
			n = 0;
		}
		break;
	}
	case 2:
	{
		for (m; m < 5; m++)   //阶段二：w*i在 sqrt(3)/3 *a——2 * sqrt(3)/3 *a
		{
			for (n; n < 8; n++)
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
				tripoint.push_back(trans(F)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(D));
				tripoint.push_back(trans(D)); tripoint.push_back(trans(A));
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
		for (m; m < 6; m++)  //阶段三：和阶段一的区别只是Y坐标全变成负
		{
			for (n; n < 8; n++)
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
				tripoint.push_back(trans(F)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(D));
				tripoint.push_back(trans(D)); tripoint.push_back(trans(A));
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
		for (m; m < 6; m++)  //阶段四：和阶段二一样
		{
			for (n; n < 8; n++)
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
				tripoint.push_back(trans(F)); tripoint.push_back(trans(C));
				tripoint.push_back(trans(C)); tripoint.push_back(trans(A));
				tripoint.push_back(trans(A)); tripoint.push_back(trans(D));
				tripoint.push_back(trans(D)); tripoint.push_back(trans(A));
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

point mindistance(point a, point b, point c)
{
	int dis1 = pow((a.x - b.x), 2) + pow((a.y - b.y), 2);
	int dis2 = pow((a.x - c.x), 2) + pow((a.y - c.y), 2);
	if (dis1 < dis2)return b;
	else return c;
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



int main(int argc, char** argv)
{

	tripoint.clear(); //把上一层的清理掉
	interpoint.clear(); 
	triangle(0);
	for (int m=0;m< tripoint.size();m+=2)  //网格的每一条线段，首尾相连线段端点重复，所以跳两个
	{
		interpoint.clear(); //上一条线段的交点清理掉
		for (int n = 0; n < polypoint.size(); n++) //多边形的每一条线段，首尾相连线段端点不重复，所以跳一个
		{
			inter(tripoint[m], tripoint[m + 1], polypoint[n], polypoint[(n + 1 + polypoint.size())% polypoint.size()]);
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
				new_tripoint.push_back(tripoint[m+1]);
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
				new_tripoint.push_back(tripoint[m+1]);
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
				new_tripoint.push_back(mindistance(tripoint[m], interpoint[0],interpoint[1]));
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




	init(argc, argv);
	glutDisplayFunc(myDisplay);
	glutMainLoop();


	return 0;

}
