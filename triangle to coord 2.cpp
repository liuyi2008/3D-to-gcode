#include <gl/glut.h>
#include<math.h>
#include <stdlib.h>//exit
#include <vector>
using std::vector;

struct point {

	double x, y;

}p1;
vector<point>tripoint; 

point trans(point p)  //坐标缩放平移  这里似乎不应该缩放，只需要平移
{
	int a = 1;  //放缩 
	int b = 0;  //左右平移
	int c = 0;  //上下平移
	p = { a * p.x + b , a * p.y + c };
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
	while (i-1 > 4* sqrt(3)*0.333*a/width)  
	{
		i = i - int(4 * sqrt(3)*0.333*a / width);
	}
	int condition;
	if (i-1< sqrt(3)*0.333*a / width)
	{
		condition = 1;
	}
	else if (i-1 < 2 * sqrt(3)*0.333*a / width )
	{
		condition = 2;
	}
	else if (i-1 < 3 * sqrt(3)*0.333*a / width )
	{
		condition = 3;
	}
	else if (i-1 < 4 * sqrt(3)*0.333*a / width )
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
		for (n; n < 5; n++)
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
		for (n; n < 5; n++)
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
		i = i - 2 * sqrt(3)*0.333*a / width ;
	for (m; m < 5; m++)  //阶段三：和阶段一的区别只是Y坐标全变成负
	{
		for (n; n < 5; n++)
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
		i = i - 2 * sqrt(3)*0.333*a / width ;
	for (m; m < 5; m++)  //阶段四：和阶段二一样
	{
		for (n; n < 5; n++)
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

	for (int i = 0; i < tripoint.size(); i++)
	{

		glColor3f(1.0f, 0.0f, 0.0f); //线是红色的
		glVertex2f(tripoint[i].x, tripoint[i].y);
		//printf("i=%d,x=%f,y=%f\n",i , tripoint[i].x, tripoint[i].y);		
		
	}
	glEnd();
	glFlush();
}

int main(int argc, char** argv)
{
	
	int i = 100;
	triangle(i);
	
	//for (int i = 0; i < tripoint.size(); i++)
	//{

	//	//glColor3f(1.0f, 0.0f, 0.0f); //线是红色的
	//	//glVertex2i(tripoint[i].x, tripoint[i].y);
	//	printf("x=%f,y=%f\n" , tripoint[i].x, tripoint[i].y);		

	//}

	init(argc, argv);
	glutDisplayFunc(myDisplay);
	glutMainLoop();
	return 0;
}