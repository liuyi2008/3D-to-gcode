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

vector<point>polypoint = { {250,50},{550,150},{550,400},{250,250},{100,350},{100,100},{120,30} }; //���������꼯
vector<point>tripoint;  //�仯�� �߶����� �ɶ�
vector<point>new_tripoint;  //�±仯�� �߶����� �ɶ�
vector<point>interpoint; //��������

//extern const string file_1;  //slice
//vector<vector<point> >model; //slice



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
		bool cond1 = (polypoint[i].x <= x0) && (x0 < polypoint[(i + 1 + n) % n].x);
		bool cond2 = (polypoint[(i + 1 + n) % n].x <= x0) && (x0 < polypoint[i].x);
		bool above = (y0 < slope * (x0 - polypoint[i].x) + polypoint[i].y);
		if ((cond1 || cond2) && above) crossings++;
	}
	return (crossings % 2 != 0);    //�� �� ֵ:  0:�� 1:��
}

point trans(point p)  //����ƽ��  �����ƺ���Ӧ�����ţ�ֻ��Ҫƽ��
{
	int b = -250;  //����ƽ��
	int c = 0;  //����ƽ��
	p = {p.x + b , p.y + c };
	return p;
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
		for (m; m < 5; m++)  //�׶�һ��w*i �� 0����sqrt(3)/3 *a
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
		for (m; m < 5; m++)   //�׶ζ���w*i�� sqrt(3)/3 *a����2 * sqrt(3)/3 *a
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
		for (m; m < 6; m++)  //�׶������ͽ׶�һ������ֻ��Y����ȫ��ɸ�
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
		for (m; m < 6; m++)  //�׶��ģ��ͽ׶ζ�һ��
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

	default: printf("������");
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



int main(int argc, char** argv)
{

	tripoint.clear(); //����һ��������
	interpoint.clear(); 
	triangle(0);
	for (int m=0;m< tripoint.size();m+=2)  //�����ÿһ���߶Σ���β�����߶ζ˵��ظ�������������
	{
		interpoint.clear(); //��һ���߶εĽ��������
		for (int n = 0; n < polypoint.size(); n++) //����ε�ÿһ���߶Σ���β�����߶ζ˵㲻�ظ���������һ��
		{
			inter(tripoint[m], tripoint[m + 1], polypoint[n], polypoint[(n + 1 + polypoint.size())% polypoint.size()]);
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
				new_tripoint.push_back(tripoint[m+1]);
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
				new_tripoint.push_back(tripoint[m+1]);
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
				new_tripoint.push_back(mindistance(tripoint[m], interpoint[0],interpoint[1]));
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




	init(argc, argv);
	glutDisplayFunc(myDisplay);
	glutMainLoop();


	return 0;

}
