#include <gl/glut.h>  
#include<math.h>
#include<stdio.h>

using namespace std;

struct point { GLfloat x, y, z; };

void main()
{
	FILE* fp;

	errno_t err;     //判断此文件流是否存在 存在返回1

	err = fopen_s(&fp, "test obj.txt", "a"); //若return 1 , 则将指向这个文件的文件流给fp1


	float a = 5.0;             //三角形边长
	float thickness = 0.8;      //壳厚
	float h = 0.5*a*0.333*sqrt(3) * 4; //45°
	//float h = 2*0.333*a;  //30°
	//float h = 2 * a;  //60°
	                                //cell高 即调整角度
								   //a、h、角度,知二得三
	int m0 = 30, n0 = 30, l0 = 20;//m列数（一行几个），n行数（有几行），l有几层

	int m = 0, n = 0, l = 0; //计数器
    float x = a + thickness, y = 0.5*sqrt(3)*(a + thickness), z = h + thickness;//偏移值 

	int r = 0;//点数计数器

	for (l = 0; l < l0; l++)
	{
		for (n = 0; n < n0; n++)
		{
			for (m = 0; m < m0; m++)
			{

				point O1 = { 0.5*a + m * x + n * 0.5*x,   0.5*0.333*sqrt(3)*a + n * y,   0.00 + l * z };

				point D2 = { 0.5*a + m * x + n * 0.5*x,   0.5*sqrt(3)*a + n * y,         0.25*h + l * z };
				point E3 = { 0.00 + m * x + n * 0.5*x,    0.00 + n * y,                  0.25*h + l * z };
				point F4 = { a + m * x + n * 0.5*x,       0.00 + n * y,                  0.25*h + l * z };

				point G5 = { 0.00 + m * x + n * 0.5*x,    a*0.33*sqrt(3) + n * y,         0.5*h + l * z };
				point H6 = { 0.5*a + m * x + n * 0.5*x,   -0.5*0.333*sqrt(3)*a + n * y,   0.5*h + l * z };
				point I7 = { a + m * x + n * 0.5*x,       a*0.33*sqrt(3) + n * y,        0.5*h + l * z };

				point D8 = { 0.5*a + m * x + n * 0.5*x,   0.5*sqrt(3)*a + n * y,          0.75*h + l * z };
				point E9 = { 0.00 + m * x + n * 0.5*x,    0.00 + n * y,                  0.75*h + l * z };
				point F10 = { a + m * x + n * 0.5*x,      0.00 + n * y,                   0.75*h + l * z };

				point O11 = { 0.5*a + m * x + n * 0.5*x,  0.5*0.333*sqrt(3)*a + n * y,   h + l * z };

				fprintf(fp, "v %f %f %f\n", O1.x, O1.y, O1.z);
				fprintf(fp, "v %f %f %f\n", D2.x, D2.y, D2.z);
				fprintf(fp, "v %f %f %f\n", E3.x, E3.y, E3.z);
				fprintf(fp, "v %f %f %f\n", F4.x, F4.y, F4.z);
				fprintf(fp, "v %f %f %f\n", G5.x, G5.y, G5.z);
				fprintf(fp, "v %f %f %f\n", H6.x, H6.y, H6.z);
				fprintf(fp, "v %f %f %f\n", I7.x, I7.y, I7.z);
				fprintf(fp, "v %f %f %f\n", D8.x, D8.y, D8.z);
				fprintf(fp, "v %f %f %f\n", E9.x, E9.y, E9.z);
				fprintf(fp, "v %f %f %f\n", F10.x, F10.y, F10.z);
				fprintf(fp, "v %f %f %f\n", O11.x, O11.y, O11.z);

				fprintf(fp, "f %d %d %d %d\n", 1 + r, 2 + r, 5 + r, 3 + r);
				fprintf(fp, "f %d %d %d %d\n", 1 + r, 3 + r, 6 + r, 4 + r);
				fprintf(fp, "f %d %d %d %d\n", 1 + r, 4 + r, 7 + r, 2 + r);
				fprintf(fp, "f %d %d %d\n", 5 + r, 2 + r, 8 + r);
				fprintf(fp, "f %d %d %d\n", 5 + r, 9 + r, 3 + r);
				fprintf(fp, "f %d %d %d\n", 6 + r, 3 + r, 9 + r);
				fprintf(fp, "f %d %d %d\n", 6 + r, 10 + r, 4 + r);
				fprintf(fp, "f %d %d %d\n", 7 + r, 4 + r, 10 + r);
				fprintf(fp, "f %d %d %d\n", 7 + r, 8 + r, 2 + r);
				fprintf(fp, "f %d %d %d %d\n", 11 + r, 9 + r, 5 + r, 8 + r);
				fprintf(fp, "f %d %d %d %d\n", 11 + r, 10 + r, 6 + r, 9 + r);
				fprintf(fp, "f %d %d %d %d\n", 11 + r, 8 + r, 7 + r, 10 + r);

				r += 11;
			}
		}
	}

	for (l = 0; l < l0 - 1; l++) //cell2，转了60°，x坐标+0.5*x,y坐标+0.5*0.333*sqrt(3)*a,z坐标+0.5*z
	{
		for (n = 0; n < n0 - 1; n++)
		{
			for (m = 0; m < m0 - 1; m++)
			{

				point O1 = { 0.5*a + m * x + n * 0.5*x + 0.5*x,   0.5*0.333*sqrt(3)*a + n * y + 0.5*0.333*sqrt(3)*a,   0.00 + l * z + 0.5*z };

				point D2 = { 0.00 + m * x + n * 0.5*x + 0.5*x,    a*0.33*sqrt(3) + n * y + 0.5*0.333*sqrt(3)*a,        0.25*h + l * z + 0.5*z };
				point E3 = { 0.5*a + m * x + n * 0.5*x + 0.5*x,   -0.5*0.333*sqrt(3)*a + n * y + 0.5*0.333*sqrt(3)*a,  0.25*h + l * z + 0.5*z };
				point F4 = { a + m * x + n * 0.5*x + 0.5*x,       a*0.33*sqrt(3) + n * y + 0.5*0.333*sqrt(3)*a,        0.25*h + l * z + 0.5*z };

				point G5 = { 0.00 + m * x + n * 0.5*x + 0.5*x,    0.00 + n * y + 0.5*0.333*sqrt(3)*a,                  0.5*h + l * z + 0.5*z };
				point H6 = { a + m * x + n * 0.5*x + 0.5*x,       0.00 + n * y + 0.5*0.333*sqrt(3)*a,                  0.5*h + l * z + 0.5*z };
				point I7 = { 0.5*a + m * x + n * 0.5*x + 0.5*x,   0.5*sqrt(3)*a + n * y + 0.5*0.333*sqrt(3)*a,         0.5*h + l * z + 0.5*z };

				point D8 = { 0.00 + m * x + n * 0.5*x + 0.5*x,    a*0.33*sqrt(3) + n * y + 0.5*0.333*sqrt(3)*a,        0.75*h + l * z + 0.5*z };
				point E9 = { 0.5*a + m * x + n * 0.5*x + 0.5*x,   -0.5*0.333*sqrt(3)*a + n * y + 0.5*0.333*sqrt(3)*a,  0.75*h + l * z + 0.5*z };
				point F10 = { a + m * x + n * 0.5*x + 0.5*x,       a*0.33*sqrt(3) + n * y + 0.5*0.333*sqrt(3)*a,       0.75*h + l * z + 0.5*z };

				point O11 = { 0.5*a + m * x + n * 0.5*x + 0.5*x,  0.5*0.333*sqrt(3)*a + n * y + 0.5*0.333*sqrt(3)*a,   h + l * z + 0.5*z };

				fprintf(fp, "v %f %f %f\n", O1.x, O1.y, O1.z);
				fprintf(fp, "v %f %f %f\n", D2.x, D2.y, D2.z);
				fprintf(fp, "v %f %f %f\n", E3.x, E3.y, E3.z);
				fprintf(fp, "v %f %f %f\n", F4.x, F4.y, F4.z);
				fprintf(fp, "v %f %f %f\n", G5.x, G5.y, G5.z);
				fprintf(fp, "v %f %f %f\n", H6.x, H6.y, H6.z);
				fprintf(fp, "v %f %f %f\n", I7.x, I7.y, I7.z);
				fprintf(fp, "v %f %f %f\n", D8.x, D8.y, D8.z);
				fprintf(fp, "v %f %f %f\n", E9.x, E9.y, E9.z);
				fprintf(fp, "v %f %f %f\n", F10.x, F10.y, F10.z);
				fprintf(fp, "v %f %f %f\n", O11.x, O11.y, O11.z);

				fprintf(fp, "f %d %d %d %d\n", 1 + r, 2 + r, 5 + r, 3 + r);
				fprintf(fp, "f %d %d %d %d\n", 1 + r, 3 + r, 6 + r, 4 + r);
				fprintf(fp, "f %d %d %d %d\n", 1 + r, 4 + r, 7 + r, 2 + r);
				fprintf(fp, "f %d %d %d\n", 5 + r, 2 + r, 8 + r);
				fprintf(fp, "f %d %d %d\n", 5 + r, 9 + r, 3 + r);
				fprintf(fp, "f %d %d %d\n", 6 + r, 3 + r, 9 + r);
				fprintf(fp, "f %d %d %d\n", 6 + r, 10 + r, 4 + r);
				fprintf(fp, "f %d %d %d\n", 7 + r, 4 + r, 10 + r);
				fprintf(fp, "f %d %d %d\n", 7 + r, 8 + r, 2 + r);
				fprintf(fp, "f %d %d %d %d\n", 11 + r, 9 + r, 5 + r, 8 + r);
				fprintf(fp, "f %d %d %d %d\n", 11 + r, 10 + r, 6 + r, 9 + r);
				fprintf(fp, "f %d %d %d %d\n", 11 + r, 8 + r, 7 + r, 10 + r);

				r += 11;
			}
		}
	}


	fprintf(fp, "v %f %f %f\n", -thickness,              -thickness - 0.5*a*0.33*sqrt(3), l0*z + thickness);
	fprintf(fp, "v %f %f %f\n", (n0 - 1)*0.5*x + m0 * x, -thickness - 0.5*a*0.33*sqrt(3), l0*z + thickness);
	fprintf(fp, "v %f %f %f\n", -thickness,              y*n0 + thickness,                l0*z + thickness);
	fprintf(fp, "v %f %f %f\n", (n0 - 1)*0.5*x + m0 * x, y*n0 + thickness,                l0*z + thickness);

	fprintf(fp, "v %f %f %f\n", -thickness, y*n0 + thickness, -thickness);
	fprintf(fp, "v %f %f %f\n", (n0 - 1)*0.5*x + m0 * x, y*n0 + thickness, -thickness);
	fprintf(fp, "v %f %f %f\n", -thickness, -thickness - 0.5*a*0.33*sqrt(3), -thickness);
	fprintf(fp, "v %f %f %f\n", (n0 - 1)*0.5*x + m0 * x, -thickness - 0.5*a*0.33*sqrt(3), -thickness);

	fprintf(fp, "f %d %d %d %d\n", 1 + r, 2 + r, 4 + r, 3 + r);
	fprintf(fp, "f %d %d %d %d\n", 7 + r, 5 + r, 6 + r, 8 + r);
	fprintf(fp, "f %d %d %d %d\n", 1 + r, 3 + r, 5 + r, 7 + r);
	fprintf(fp, "f %d %d %d %d\n", 8 + r, 6 + r, 4 + r, 2 + r);
	fprintf(fp, "f %d %d %d %d\n", 7 + r, 8 + r, 2 + r, 1 + r);
	fprintf(fp, "f %d %d %d %d\n", 3 + r, 4 + r, 6 + r, 5 + r);

	fclose(fp);

}