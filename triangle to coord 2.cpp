#include <gl/glut.h>
#include<math.h>
#include <stdlib.h>//exit
#include <vector>
using std::vector;

struct point {

	double x, y;

}p1;
vector<point>tripoint; 

point trans(point p)  //��������ƽ��  �����ƺ���Ӧ�����ţ�ֻ��Ҫƽ��
{
	int a = 1;  //���� 
	int b = 0;  //����ƽ��
	int c = 0;  //����ƽ��
	p = { a * p.x + b , a * p.y + c };
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
	//************************�������i,ʹ֮ӳ��������0��4*sqrt(3)*0.333*a��*****************************
	switch (condition)
	{
	case 1:
	{
	for (m; m < 5; m++)  //�׶�һ��w*i �� 0����sqrt(3)/3 *a
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
	for (m; m < 5; m++)   //�׶ζ���w*i�� sqrt(3)/3 *a����2 * sqrt(3)/3 *a
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
	for (m; m < 5; m++)  //�׶������ͽ׶�һ������ֻ��Y����ȫ��ɸ�
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
	for (m; m < 5; m++)  //�׶��ģ��ͽ׶ζ�һ��
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

	default: printf("������");
		break;
	}

}

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

	for (int i = 0; i < tripoint.size(); i++)
	{

		glColor3f(1.0f, 0.0f, 0.0f); //���Ǻ�ɫ��
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

	//	//glColor3f(1.0f, 0.0f, 0.0f); //���Ǻ�ɫ��
	//	//glVertex2i(tripoint[i].x, tripoint[i].y);
	//	printf("x=%f,y=%f\n" , tripoint[i].x, tripoint[i].y);		

	//}

	init(argc, argv);
	glutDisplayFunc(myDisplay);
	glutMainLoop();
	return 0;
}