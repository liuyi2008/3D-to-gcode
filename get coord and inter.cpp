#include<stdio.h>
#include<iostream>
#include <vector>
#include<math.h>
#include<windows.h>
#include<GL/glut.h>
using std::vector;
using namespace std;
struct point {

	float x, y;

}p1;

typedef struct XET
{
	float x;
	float dx, ymax;
	XET* next;
} AET, NET;

const int POINTNUM = 7;      //����ε���.
vector<point>polypoint = { {250,50},{550,150},{550,400},{250,250},{100,350},{100,100},{120,30} }; //���������꼯
vector<vector<point> >coord1;   //�仯�ߵ����꼯, һԪ����һ��vector��װһ��������    //vector<point>coord1;  //�仯�ߵ����꼯
vector<point>coord;    //��coord1���
vector<point>coord2;  //��������꼯
vector<vector<point> >coord3; //ȡ��������꼯    
//vector<point>coord4;  

void getcoord(point A, point B)  //�������ĵ����coord��
{
	float xmin, xmax, y, ymin, ymax;
	float x1, x2, y1, y2;
	x1 = A.x;
	y1 = A.y;
	x2 = B.x;
	y2 = B.y;

	float a = (y2 - y1) / (x2 - x1);
	float b = y1 - x1 * a;

	if (A.x > B.x)
	{
		xmin = x2;
		ymin = x2;
		xmax = x1;
		ymax = y1;
	}
	else
	{
		xmin = x1;
		ymin = y1;
		xmax = x2;
		ymax = y2;

	}

	for (; xmin < xmax; xmin++)
	{
		y = a * xmin + b;
		p1 = { roundf(xmin) , roundf(y) };
		coord.push_back(p1);

	}
	p1 = { roundf(xmax) , roundf(ymax) };
	coord.push_back(p1);

}

/*
bool isNotMember(point a, int j)   //�ǲ���coord1[j]�ĳ�Ա
{
	bool f = false;
	for (int i = 0; i < coord1[j].size(); i++)
		if (coord1[j][i].x == a.x && coord1[j][i].y == a.y)
		{
			f = true;
			break;
		}
	return f;
}
void inter() //��1��2���еķŵ�3
{
	
	for (int i = 0; i < coord2.size(); i++)
		for (int j = 0; j < ; i++)
		{
			if (isNotMember(coord2[i],j))
			{
				coord3.push_back(coord2[i]);
			}
		}
		
}
*/

bool isNotMember(point a)   //�ǲ���coord2�ĳ�Ա
{
	bool f = false;
	for (int i = 0; i < coord2.size(); i++)
		if (coord2[i].x == a.x && coord2[i].y == a.y)
		{
			f = true;
			break;
		}
	return f;
}
void inter() //��1��2���еķŵ�3
{

	for (int i = 0; i < coord1.size(); i++)
	{
		for (int j = 0; j < coord1[i].size(); j++)
		{
			
			if (isNotMember(coord1[i][j]))
			{
				p1 = { coord1[i][j].x , coord1[i][j].y };
				printf("x=%f,y=%f,i=%d,j=%d\n", coord1[i][j].x , coord1[i][j].y,i,j);
				coord.push_back(p1);
			}
		}

		coord3.push_back(coord);
		coord.clear();

	}
}

void PolyScan()   //��������polypoint�õ��������꼯coord2���˺����е�i�����i�л���y����
{
	/******������ߵ��y����(ɨ�赽�˽���)****************************************/
	int MaxY = 0;
	int i;
	for (i = 0; i < POINTNUM; i++)
		if (polypoint[i].y > MaxY)
			MaxY = polypoint[i].y;

	/*******��ʼ��AET��***********************************************************/
	AET *pAET = new AET;
	pAET->next = NULL;

	/******��ʼ��NET��************************************************************/
	NET *pNET[1024];
	for (i = 0; i <= MaxY; i++)
	{
		pNET[i] = new NET;
		pNET[i]->next = NULL;
	}
	/******ɨ�貢����NET��*********************************************************/
	for (i = 0; i <= MaxY; i++)
	{
		for (int j = 0; j < POINTNUM; j++)
			if (polypoint[j].y == i)
			{
				//һ�����ǰ���һ�����γ�һ���߶Σ�������ĵ�Ҳ�γ��߶�
				if (polypoint[(j - 1 + POINTNUM) % POINTNUM].y > polypoint[j].y)
				{
					NET *p = new NET;
					p->x = polypoint[j].x;
					p->ymax = polypoint[(j - 1 + POINTNUM) % POINTNUM].y;
					p->dx = (polypoint[(j - 1 + POINTNUM) % POINTNUM].x - polypoint[j].x) / (polypoint[(j - 1 + POINTNUM) % POINTNUM].y - polypoint[j].y);
					p->next = pNET[i]->next;
					pNET[i]->next = p;
				}
				if (polypoint[(j + 1 + POINTNUM) % POINTNUM].y > polypoint[j].y)
				{
					NET *p = new NET;
					p->x = polypoint[j].x;
					p->ymax = polypoint[(j + 1 + POINTNUM) % POINTNUM].y;
					p->dx = (polypoint[(j + 1 + POINTNUM) % POINTNUM].x - polypoint[j].x) / (polypoint[(j + 1 + POINTNUM) % POINTNUM].y - polypoint[j].y);
					p->next = pNET[i]->next;
					pNET[i]->next = p;
				}
			}
	}
	/******���������»��Ա߱�AET*****************************************************/
	for (i = 0; i <= MaxY; i++)
	{
		//�����µĽ���x,����AET
		NET *p = pAET->next;
		while (p)
		{
			p->x = p->x + p->dx;
			p = p->next;
		}
		//���º���AET������*************************************************************/
		//�ϱ�����,���ٿ��ٿռ�
		AET *tq = pAET;
		p = pAET->next;
		tq->next = NULL;
		while (p)
		{
			while (tq->next && p->x >= tq->next->x)
				tq = tq->next;
			NET *s = p->next;
			p->next = tq->next;
			tq->next = p;
			p = s;
			tq = pAET;
		}
		//(�Ľ��㷨)�ȴ�AET����ɾ��ymax==i�Ľ��****************************************/
		AET *q = pAET;
		p = q->next;
		while (p)
		{
			if (p->ymax == i)
			{
				q->next = p->next;
				delete p;
				p = q->next;
			}
			else
			{
				q = q->next;
				p = q->next;
			}
		}
		//��NET�е��µ����AET,���ò��뷨��Xֵ��������**********************************/
		p = pNET[i]->next;
		q = pAET;
		while (p)
		{
			while (q->next && p->x >= q->next->x)
				q = q->next;
			NET *s = p->next;
			p->next = q->next;
			q->next = p;
			p = s;
			q = pAET;
		}
		/******��������ɫ***************************************************************/

		p = pAET->next;
		while (p && p->next)
		{
			for (float j = p->x; j <= p->next->x; j++)
			{
				//glVertex2i(static_cast<int>(j), i);
				//printf("coord2:x:%f,y=%f\n", roundf(j), float(i));
				p1 = { roundf(j) , float(i) };
				coord2.push_back(p1);
			}
			p = p->next->next;//���Ƕ˵����
		}
	}
}

point trans(point p)  //��������ƽ��
{
	int a = 100; //����
	int b = 0; //ƽ��
 	p = { a * p.x + b , a * p.y + b };
	return p;
}

void triangle()
{
	float Thickness = 0.2; // һ��ĺ��

	float a = 1.0f;//��λ����

	float width = 2; //��ͷ�Ŀ��

	int i = 0; // ��i�㣬��i��wһ�������A��B��C��λ��

	int m = 0;  //��
	int n = 0;  //��

	for (m; m < 5; m++)
	{
		for (n; n < 5; n++)
		{
			point A = { m*0.5*a + n * a + 0.5*a ,                 m*0.5*sqrt(3)*a + sqrt(3) * 0.333 * (3 * 0.5 * a - 2 * sqrt(3)*width*i) };

			point B = { m*0.5*a + n * a + sqrt(3)*width*i ,       m*0.5*sqrt(3)*a + width * i };

			point C = { m*0.5*a + n * a + a - sqrt(3)*width*i ,   m*0.5*sqrt(3)*a + width * i };

			point D = { m*0.5*a + n * a + 0.5*a ,                 m*0.5*sqrt(3)*a + 0.5*sqrt(3)*a };

			point E = { m * 0.5 * a + n * a + 0 ,                 m*0.5*sqrt(3)*a + 0 };

			point F = { m*0.5*a + n * a + a ,                     m*0.5*sqrt(3)*a + 0 };

			getcoord(trans(E), trans(B));
			getcoord(trans(B), trans(C));
			getcoord(trans(C), trans(F));
			getcoord(trans(F), trans(C));
			getcoord(trans(C), trans(A));
			getcoord(trans(A), trans(D));
			getcoord(trans(D), trans(A));
			getcoord(trans(A), trans(B));
			getcoord(trans(B), trans(E));

			coord1.push_back(coord);
			coord.clear();
		}
		n = 0;
	}
}


//**********************************����������*******************
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
	glBegin(GL_POINTS);

	for (int i = 0; i < coord1.size(); i++)
	{
		for (int j = 0; j < coord1[i].size(); j++)
		{
			glColor3f(1.0f, 0.0f, 0.0f); //���Ǻ�ɫ��
			glVertex2i(coord1[i][j].x, coord1[i][j].y);
			//printf("x=%d,y=%d\n", coord1[i][j].x, coord1[i][j].y);
		}
	}

	for (int i = 0; i < coord2.size(); i++)
	{
		glColor3f(0.0f, 0.0f, 1.0f); //���������ɫ��
		glVertex2i(coord2[i].x, coord2[i].y);
	}

	for (int i = 0; i < coord3.size(); i++)
	{
		for (int j = 0; j < coord3[i].size(); j++)
		{
			glColor3f(0.0f, 1.0f, 0.0f); //������ɫ��
			glVertex2i(coord3[i][j].x, coord3[i][j].y);
			//printf("x=%d,y=%d\n", coord3[i][j].x, coord3[i][j].y);
		}
	}

	glEnd();
	glFlush();
}

int main(int argc, char** argv)
{

	triangle();
	PolyScan();
	//inter();
	printf("%d", coord2.size());

	init(argc, argv);
	glutDisplayFunc(myDisplay);        //ͼ�εĶ��崫�ݸ���window.
	glutMainLoop();
	return 0;
}

//**********************************����������*******************
/*
int main() 
{
	point p1 = { 123 , 452 };
	point p2 = { 457 , -11 };

	getcoord1(p1, p2);  //�������㣬���������������굽coord1
	for (int i = 0; i < coord1.size(); i++)
	{
		printf("coord1:x=%f,y=%f\n", coord1[i].x, coord1[i].y);
	}

	PolyScan();//����������polypoint������������꼯coord2

	inter();//ȡ������1��2���еķŵ�3

	for (int i = 0; i < coord3.size(); i++)
	{
		printf("coord3:x=%f,y=%f\n", coord3[i].x, coord3[i].y);

	}
	return 0;
}
*/


//�Ľ����룺�ܲ����ñ仯��������Σ�

