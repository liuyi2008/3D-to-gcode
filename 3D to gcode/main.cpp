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

vector<point>polypoint; //һ�������������
vector<point>tripoint;  //һ��ı仯�� �߶����� �ɶ�
vector<point>new_tripoint;  //�µ�һ��ı仯�� �߶����� �ɶ�
vector<point>interpoint; //��������

extern const string file_1;  //slice
extern const string file_2; 

vector<vector<point> >model; //slice  ����ģ�͵ģ��ֲ�����������꣬һά��������ά�˲�������ߵ�����
vector<vector<point> >modelfill(model); //��ɵ�model����ߣ�һά��������ά�˲�ı仯�� �߶�����

//**************************************************************�ж��߶����߶εĹ�ϵ����********************************************************
//���㽻��˻�(P1-P0)x(P2-P0)
double xmult(point p1, point p2, point p0) {
	return (p1.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p1.y - p0.y);
}

//�е��Ƿ����߶���,�����˵�
int dot_online_in(point p, point l1, point l2) {
	return zero(xmult(p, l1, l2)) && (l1.x - p.x)*(l2.x - p.x) < eps && (l1.y - p.y)*(l2.y - p.y) < eps;
}

//���������߶�ͬ��,�����߶��Ϸ���0
int same_side(point p1, point p2, point l1, point l2) {
	return xmult(l1, p1, l2)*xmult(l1, p2, l2) > eps;
}

//����ֱ��ƽ��
int parallel(point u1, point u2, point v1, point v2) {
	return zero((u1.x - u2.x)*(v1.y - v2.y) - (v1.x - v2.x)*(u1.y - u2.y));
}

//�����㹲��
int dots_inline(point p1, point p2, point p3) {
	return zero(xmult(p1, p2, p3));
}

//�����߶��ཻ,�����˵�Ͳ����غ�
int intersect_in(point u1, point u2, point v1, point v2) {
	if (!dots_inline(u1, u2, v1) || !dots_inline(u1, u2, v2))
		return !same_side(u1, u2, v1, v2) && !same_side(v1, v2, u1, u2);
	return dot_online_in(u1, v1, v2) || dot_online_in(u2, v1, v2) || dot_online_in(v1, u1, u2) || dot_online_in(v2, u1, u2);
}

//�������߶ν���,�����߶��Ƿ��ཻ(ͬʱ����Ҫ�ж��Ƿ�ƽ��!)
point intersection(point u1, point u2, point v1, point v2) {
	point ret = u1;
	double t = ((u1.x - v1.x)*(v1.y - v2.y) - (u1.y - v1.y)*(v1.x - v2.x))
		/ ((u1.x - u2.x)*(v1.y - v2.y) - (u1.y - u2.y)*(v1.x - v2.x));
	ret.x += (u2.x - u1.x)*t;
	ret.y += (u2.y - u1.y)*t;
	return ret;
}

void inter(point u1, point u2, point v1, point v2)  //���������н���潻�㣬û�����޲���
{
	point ans;
	if (parallel(u1, u2, v1, v2) || !intersect_in(u1, u2, v1, v2))
	{

	}
	else {
		ans = intersection(u1, u2, v1, v2);
		interpoint.push_back(ans);
		//printf("����Ϊ:(%lf,%lf)", ans.x, ans.y);
	}

}

//****************************************************************�ж��߶����߶εĹ�ϵ����************************************************************

int InOrOutPolygon(point a)  //�жϵ��Ƿ��ڶ������
{
	double x0 = a.x;
	double y0 = a.y;
	int crossings = 0;
	int n = polypoint.size();
	for (int i = 0; i < n; i++)
	{
		// ��������x֮�� ���Ե㴹ֱy������������
		double slope = (polypoint[(i + 1 + n) % n].y - polypoint[i].y) / (polypoint[(i + 1 + n) % n].x - polypoint[i].x);
		boolean cond1 = (polypoint[i].x <= x0) && (x0 < polypoint[(i + 1 + n) % n].x);
		boolean cond2 = (polypoint[(i + 1 + n) % n].x <= x0) && (x0 < polypoint[i].x);
		boolean above = (y0 < slope * (x0 - polypoint[i].x) + polypoint[i].y);
		if ((cond1 || cond2) && above) crossings++;
	}
	return (crossings % 2 != 0);    //�� �� ֵ:  0:�� 1:��
}

point trans(point p)  //����ƽ��  �����ƺ���Ӧ�����ţ�ֻ��Ҫƽ��
{
	int b = -250;  //����ƽ��
	int c = 0;  //����ƽ��
	p = { p.x + b , p.y + c };
	return p;
}

point trans2(point p)//����������
{
	int a = 1;  //����
	int b = -250;  //����ƽ��
	int c = 0;  //����ƽ��
	p = { a * p.x + b , a * p.y + c };
	return p;
}

point mindistance(point a, point b, point c) //��b �� c�з�����a ����ĵ�
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

void triangle(int i)  //��i�㣬��1��ʼ����i��wһ�������A��B��C��λ��
{
	float Thickness = 0.2; // һ��ĺ��

	float a = 100.0f;//��λ���� ��������������

	float width = 0.2; //��ͷ�Ŀ��

	int m = 0;  //��
	int n = 0;  //��
	//************************�������i,ʹ֮ӳ��������0��4*sqrt(3)*0.333*a��*****************************
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
	//************************�������i,ʹ֮ӳ��������0��4*sqrt(3)*0.333*a��*****************************
	switch (condition)
	{
	case 1:
	{
		for (m; m < 10; m++)  //�׶�һ��w*i �� 0����sqrt(3)/3 *a
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
		for (m; m < 10; m++)   //�׶ζ���w*i�� sqrt(3)/3 *a����2 * sqrt(3)/3 *a
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
		for (m; m < 10; m++)  //�׶������ͽ׶�һ������ֻ��Y����ȫ��ɸ�
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
		for (m; m < 10; m++)  //�׶��ģ��ͽ׶ζ�һ��
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

	default: printf("������");
		break;
	}

}

//****************************************************************������������************************************************************

void init(int argc, char** argv)
{
	glutInit(&argc, argv);  //I��ʼ�� GLUT.
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);  //������ʾģʽ�����������ʹ��RGBģ��
	glutInitWindowPosition(50, 100);  //���ô��ڵĶ��������λ��
	glutInitWindowSize(400, 300);  //���ô��ڵĸ߶ȺͿ��
	glutCreateWindow("Scan Program");

	glClearColor(1.0, 1.0, 1.0, 0); //���ڱ�����ɫΪ��ɫ
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

		glColor3f(1.0f, 0.0f, 0.0f); //���Ǻ�ɫ��
		glVertex2f(new_tripoint[i].x, new_tripoint[i].y);
		//printf("i=%d,x=%f,y=%f\n",i , tripoint[i].x, tripoint[i].y);		

	}

	glEnd();
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < polypoint.size(); i++)
	{

		glColor3f(0.0f, 0.0f, 1.0f); //���Ǻ�ɫ��
		glVertex2f(polypoint[i].x, polypoint[i].y);
		//printf("i=%d,x=%f,y=%f\n",i , tripoint[i].x, tripoint[i].y);		

	}
	glEnd();
	glFlush();
}

//****************************************************************������������************************************************************

void main(int argc, char** argv)
{
	printf("slicing...\n");
	readfile(file_2);
	BoundingBox();
	findIntersect();  //ǰ����������Ƭ������model

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
		tripoint.clear(); //����һ��ı仯�������
		interpoint.clear(); //����һ��Ľ��������
		polypoint.clear();//����һ��������������
		triangle(i); //������i��ı仯�ߴ���tripoint
		for (int i = 0; i < tripoint.size(); i++) //�Ա仯������
		{

			p1 = { tripoint[i].x*0.2,tripoint[i].y*0.2};
			tripoint[i] = p1;

		}
		polypoint = model[i];//�ѵ�i��������߸�polypoint

		for (int m = 0; m < tripoint.size(); m += 2)  //�����ÿһ���߶Σ���β�����߶ζ˵��ظ�������������
		{
			interpoint.clear(); //��һ���߶εĽ��������
			for (int n = 0; n < polypoint.size(); n++) //����ε�ÿһ���߶Σ���β�����߶ζ˵㲻�ظ���������һ��
			{
				inter(tripoint[m], tripoint[m + 1], polypoint[n], polypoint[(n + 1 + polypoint.size()) % polypoint.size()]);
			}
			int a = InOrOutPolygon(tripoint[m]);
			int b = InOrOutPolygon(tripoint[m + 1]);  //�ж��ǲ����ڵ�,0:�� 1:��

			if (interpoint.size() == 0)   // �޽���
			{
				if ((a == 0) && (b == 0)) //���ڵ㣬��
				{
					continue;
				}
				if ((a == 1) || (b == 1)) //�����ڵ㣬ȡ
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(tripoint[m + 1]);
					continue;
				}
			}

			if (interpoint.size() == 1)   //��һ����
			{
				if ((a == 0) && (b == 0))   //��һ���������ڵ㣬�� �����У�
				{
					continue;
				}
				if ((a == 1) && (b == 0))//��һ��������һ�ڵ㣬ȡ�ڵ㵽���� ���ཻ��
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(interpoint[0]);
					continue;
				}
				if ((a == 0) && (b == 1))//��һ��������һ�ڵ㣬ȡ�ڵ㵽���� ���ཻ��
				{
					new_tripoint.push_back(interpoint[0]);
					new_tripoint.push_back(tripoint[m + 1]);
					continue;
				}
				if ((a == 1) && (b == 1)) //��һ�����������ڵ㣬ȡ ���ڲ�������
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(tripoint[m + 1]);
					continue;
				}

			}

			if (interpoint.size() == 2)   //��������
			{
				if ((a == 0) && (b == 0))   //�������������ڵ㣬ȡ������
				{
					new_tripoint.push_back(interpoint[0]);
					new_tripoint.push_back(interpoint[1]);
					continue;
				}
				if ((a == 1) && (b == 0))//����������һ���ڵ� ����������ע��˳�� 
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(mindistance(tripoint[m], interpoint[0], interpoint[1]));
					continue;
				}
				if ((a == 0) && (b == 1))//����������һ���ڵ� ����������ע��˳��
				{
					new_tripoint.push_back(tripoint[m + 1]);
					new_tripoint.push_back(mindistance(tripoint[m + 1], interpoint[0], interpoint[1]));
					continue;
				}
				if ((a == 1) && (b == 1))//���������������ڵ� ����������ע��˳�� 
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

	printf("path planning complete\n"); //modelfill���


	/*init(argc, argv);
	glutDisplayFunc(myDisplay);
	glutMainLoop();*/

	//*******************************************�����Ǳ�дGcode*******************************************

	printf("Gcode writing...\n");

	FILE* fp;

	errno_t err;     //�жϴ��ļ����Ƿ���� ���ڷ���1

	err = fopen_s(&fp, "test gcode.txt", "a"); //��return 1 , ��ָ������ļ����ļ�����fp1

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
	for (int j = 0; j < 18; j++) //����
	{
		fprintf(fp, "G0 F300 X%.3f Y%.3f Z%.3f\n", j * 4 + 30.00, 30.00,  0.500);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 30.00, 100.00, E += distance({ j * 4 + 30.00, 30.00 }, { j * 4 + 30.00, 100.00 })*0.3);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 32.00, 100.00, E += distance({ j * 4 + 30.00, 100.00 }, { j * 4 + 32.00, 100.00 })*0.3);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 32.00, 30.00,  E += distance({ j * 4 + 32.00, 100.00 }, { j * 4 + 32.00, 30.00 })*0.3);
		fprintf(fp, "G1 F300 X%.3f Y%.3f E%.3f\n", j * 4 + 34.00, 30.00,  E += distance({ j * 4 + 32.00, 30.00 }, { j * 4 + 34.00, 30.00 })*0.3);

	}

	//for (int i = 0; i < 1; i++) //��һ���Ĵ�һЩ���൱�ڵ���
	//{
	//	fprintf(fp, ";LAYER:%d\n", i);

	//	fprintf(fp, "G0 F1000 X%.3f Y%.3f Z%.3f\n", 0.00, 0.00, 0.400 + i * 0.200);  //ԭ��
	//	fprintf(fp, "G1 F1000 X%.3f Y%.3f E%.3f\n", model[i][0].x, model[i][0].y, E += distance({0.0}, model[i][0])*0.6);
	//	fprintf(fp, ";TYPE:OUTLINE\n");
	//	for (int j = 1; j < model[i].size(); j++)
	//	{
	//		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][j].x, model[i][j].y, E += distance(model[i][j - 1], model[i][j])*0.6);
	//	}
	//	fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][0].x, model[i][0].y, E += distance(model[i][model[i].size() - 1], model[i][0])*0.6); //��һ��Ȧ��Ҫ��ԭ��
	//	fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][model[i].size() - 1], model[i][0])*0.01);

	//	fprintf(fp, ";TYPE:FILL\n");
	//	for (int k = 1; k < modelfill[i].size(); k += 2)
	//	{
	//		fprintf(fp, "G0 F2600 X%.3f Y%.3f\n", modelfill[i][k - 1].x, modelfill[i][k - 1].y);
	//		fprintf(fp, "G1 F1000 X%.3f Y%.3f E%.5f\n", modelfill[i][k].x, modelfill[i][k].y, E += distance(modelfill[i][k - 1], modelfill[i][k])*0.1);
	//	}
	//}

	for (int i = 1; i < model.size(); i++) //ÿһ�� i�������
	{
		fprintf(fp, ";LAYER:%d\n", i);	
		fprintf(fp, "G0 F2000 X%.3f Y%.3f Z%.3f\n", model[i][0].x, model[i][0].y, 0.500 + i * 0.300);
		fprintf(fp, ";TYPE:OUTLINE\n");
		fprintf(fp, "G1 F1000 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][0], model[i][1])*0.1);
		for (int j = 2; j < model[i].size(); j++)
		{
			fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n",model[i][j].x, model[i][j].y, E += distance(model[i][j-1], model[i][j])*0.1);
		}
		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][0].x, model[i][0].y, E += distance(model[i][model[i].size()-1], model[i][0])*0.1); //��һ��Ȧ��Ҫ��ԭ��

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
