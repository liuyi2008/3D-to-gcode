#include <gl/glut.h>
#include<math.h>
#include <stdlib.h>//exit
#include <vector>


float Thickness = 0.2;

float a = 1.0f;//单位距离

float width = 0.2;
/*
struct A
{
	float x = 0.5*a;
	float y = sqrt(3) / 3 * (3 / 2 * a - 2 * sqrt(3)*width*i);
}A;

struct B
{
	float x = sqrt(3)*width*i;
	float y = width * i;
}B;

struct C
{
	float x = a - sqrt(3)*width*i;
	float y = width * i;
}C;
*/
struct A
{
	float x = 0.5*a;

}A;
float A_y(int i) 
{
	float y = sqrt(3) / 3 * (3 / 2 * a - 2 * sqrt(3)*width*i);
	return y;
}

float B_x(int i) 
{
	float x = sqrt(3)*width*i;
	return x;
}

float B_y(int i)
{
	float y = width * i;
	return y;
}

float C_x(int i)
{
	float x = a - sqrt(3)*width*i;
	return x;
}

float C_y(int i)
{
	float y = width * i;
	return y;
}

struct D
{
	float x = 0.5*a;
	float y = 0.5*sqrt(3)*a;
}D;

struct E
{
	float x = 0;
	float y = 0;
}E;

struct F
{
	float x = a;
	float y = 0;
}F;

float trans_x(float old)
{
	float new_ncoordinate;
	new_ncoordinate = 3 * old + 85;
	return new_ncoordinate;

}

float trans_y(float old)
{
	float new_ncoordinate;
	new_ncoordinate = -(3 * old + 40);
		return new_ncoordinate;

}

int main()
{

	FILE* fp;

	errno_t err;     //判断此文件流是否存在 存在返回1

	err = fopen_s(&fp, "test gcode.txt", "a"); //若return 1 , 则将指向这个文件的文件流给fp1

	int i = 3; //层数

	int j = 0, l = 0;//行，列数

for(i;i<100;i++)
{
	fprintf(fp, "M73 P%d\n", i);

	
	int j;
	switch (j = i % 2)
	{
	
	case 0: 
	{
		fprintf(fp, "G1 X%f Y%f Z%f F1500;travel move\n", trans_x(E.x), trans_y(E.y), i*Thickness); //回原点

		for (int t = 0; t < 10; t++)                        //第一层
		{

			fprintf(fp, "G1 X%f Y%f Z%f F1500;Inset\n", trans_x(E.x), trans_y(4 * t * width + E.y), i*Thickness);

			fprintf(fp, "G1 X%f Y%f Z%f F1500;Inset\n", trans_x(10 * a + E.x), trans_y(4 * t * width + E.y), i*Thickness);

			fprintf(fp, "G1 X%f Y%f Z%f F1500;Inset\n", trans_x(10 * a + E.x), trans_y(4 * t * width + width * 2 + E.y), i*Thickness);

			fprintf(fp, "G1 X%f Y%f Z%f F1500;Inset\n", trans_x(E.x), trans_y(4 * t * width + width * 2 + E.y), i*Thickness);

		}
	}

		//i = 2;
	case 1:
	{
		fprintf(fp, "G1 X%f Y%f Z%f F1500;travel move\n", trans_x(E.x), trans_y(E.y), i*Thickness); //回原点

		for (int t = 0; t < 10; t++)   //第二层
		{

			fprintf(fp, "G1 X%f Y%f Z%f F1500;Inset\n", trans_x(8 * t * width + E.x), trans_y(E.y), i*Thickness);

			fprintf(fp, "G1 X%f Y%f Z%f F1500;Inset\n", trans_x(8 * t * width + E.x), trans_y(5 * a + E.y), i*Thickness);

			fprintf(fp, "G1 X%f Y%f Z%f F1500;Inset\n", trans_x(8 * t * width + 2 * width + E.x), trans_y(5 * a + E.y), i*Thickness);

			fprintf(fp, "G1 X%f Y%f Z%f F1500;Inset\n", trans_x(8 * t * width + 2 * width + E.x), trans_y(E.y), i*Thickness);

		}
	}

	}

	//头两层打底座
}
//	i = 3;

	for (i; i <= 30; i++) //30个layer
	{
		
		fprintf(fp, "M73 P%d\n", i); 		
		fprintf(fp, "G1 X%f Y%f Z%f F1500;travel move\n", trans_x(j*0.5*a + l * a + E.x), trans_y(j*0.5*sqrt(3)*a + E.y), i*Thickness);
		//fprintf(fp, "M126 T0; fan on\n");

		for (j; j < 3; j++)        //3列
		{
			
			for (l;  l< 3; l++)         //3行
			{
				fprintf(fp, "G1 X%f Y%f Z%f F9000;travel move\n",trans_x(j*0.5*a+l*a+E.x),trans_y(j*0.5*sqrt(3)*a+E.y),i*Thickness);
				printf("x:%f,y:%f,层数i=%d\n", trans_x(j*0.5*a+l*a + E.x), trans_y(j*0.5*sqrt(3)*a + E.y), i);

				fprintf(fp, "G1 X%f Y%f Z%f F1800;infill\n", trans_x(j*0.5*a + l * a + E.x), trans_y(j*0.5*sqrt(3)*a + E.y), i*Thickness);
				printf("x:%f,y:%f,层数i=%d\n", trans_x(j*0.5*a + l * a + E.x), trans_y(j*0.5*sqrt(3)*a + E.y), i);;

				fprintf(fp, "G1 X%f Y%f Z%f F1800;infill\n", trans_x(j*0.5*a + l * a + B_x(i)), trans_y(j*0.5*sqrt(3)*a + B_y(i)), i*Thickness);
				printf("x:%f,y:%f,层数i=%d\n", trans_x(j*0.5*a + l * a + B_x(i)), trans_y(j*0.5*sqrt(3)*a + B_y(i)), i);

				fprintf(fp, "G1 X%f Y%f Z%f F1800;infill\n", trans_x(j*0.5*a + l * a + C_x(i)), trans_y(j*0.5*sqrt(3)*a + C_y(i)), i*Thickness);
				printf("x:%f,y:%f,层数i=%d\n", trans_x(j*0.5*a + l * a + C_x(i)), trans_y(j*0.5*sqrt(3)*a + C_y(i)), i);

				fprintf(fp, "G1 X%f Y%f Z%f F1800;infill\n", trans_x(j*0.5*a + l * a + F.x), trans_y(j*0.5*sqrt(3)*a + F.y), i*Thickness);
				printf("x:%f,y:%f,层数i=%d\n", trans_x(j*0.5*a + l * a + F.x), trans_y(j*0.5*sqrt(3)*a + F.y), i);

				fprintf(fp, "G1 X%f Y%f Z%f F1800;infill\n", trans_x(j*0.5*a + l * a + C_x(i)), trans_y(j*0.5*sqrt(3)*a + C_y(i)), i*Thickness);
				printf("x:%f,y:%f,层数i=%d\n", trans_x(j*0.5*a + l * a + C_x(i)), trans_y(j*0.5*sqrt(3)*a + C_y(i)), i);

				fprintf(fp, "G1 X%f Y%f Z%f F1800;infill\n", trans_x(j*0.5*a + l * a + A.x), trans_y(j*0.5*sqrt(3)*a + A_y(i)), i*Thickness);
				printf("x:%f,y:%f,层数i=%d\n", trans_x(j*0.5*a + l * a + A.x), trans_y(j*0.5*sqrt(3)*a + A_y(i)), i);

				fprintf(fp, "G1 X%f Y%f Z%f F1800;infill\n", trans_x(j*0.5*a + l * a + D.x), trans_y(j*0.5*sqrt(3)*a + D.y), i*Thickness);
				printf("x:%f,y:%f,层数i=%d\n", trans_x(j*0.5*a + l * a + D.x), trans_y(j*0.5*sqrt(3)*a + D.y), i);

				fprintf(fp, "G1 X%f Y%f Z%f F1800;infill\n", trans_x(j*0.5*a + l * a + A.x), trans_y(j*0.5*sqrt(3)*a + A_y(i)), i*Thickness);
				printf("x:%f,y:%f,层数i=%d\n", trans_x(j*0.5*a + l * a + A.x), trans_y(j*0.5*sqrt(3)*a + A_y(i)), i);

				fprintf(fp, "G1 X%f Y%f Z%f F1800;infill\n", trans_x(j*0.5*a + l * a + B_x(i)), trans_y(j*0.5*sqrt(3)*a + B_y(i)), i*Thickness);
				printf("x:%f,y:%f,层数i=%d\n", trans_x(j*0.5*a + l * a + B_x(i)), trans_y(j*0.5*sqrt(3)*a + B_y(i)), i);
			}

			l = 0;
		}

		//fprintf(fp, "M127 T0; fan off\n");

		j = 0;

	}

	fclose(fp);

	return 0;

}


