#include <iostream>
#include <fstream>
#include<string>
using namespace std;

struct point { float x, y, z; };

int main() 
{
	ofstream File;

	float a0 = 4;         //初始边长
	float thickness = 0.5;//初始壁厚
	float length = 150;
	float width = 150;
	float high = 100;     //盒子长高宽
	
	//int m0 = 8, n0 = 8, l0 = 8;

	for (float a = a0; a < 23; a+=2)
	{
		for (float t = thickness; t < 3.5; t+=0.5)
		{
			char cc[1024];
			sprintf_s(cc, "a=%.1ft=%.1f.obj",a,t);
			File.open(cc);

			float h = 0.5*a*0.333*sqrt(3) * 4;

			int m = 0, n = 0, l = 0;
			float x = a + t, y = 0.5*sqrt(3)*(a + t), z = h + t;
			int r = 0;
			//point O1 = { 0.5*a + m * x + n * 0.5*x,   0.5*0.333*sqrt(3)*a + n * y,   0.00 + l * z };
			//sprintf_s(cc, "v %f %f %f\n", O1.x, O1.y, O1.z);
			//File << (cc);
			//sprintf_s(cc, "v %f %f %f\n", O1.x+1, O1.y+2, O1.z+3);
			//File << (cc);

			for (l = 0; h + l * z < high; l++)
			{
				for (n = 0; 0.5*sqrt(3)*a + n * y < width; n++)
				{
					for (m = 0; a + m * x + n * 0.5*x < length; m++)
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

						sprintf_s(cc, "v %f %f %f\n", O1.x, O1.y, O1.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", D2.x, D2.y, D2.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", E3.x, E3.y, E3.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", F4.x, F4.y, F4.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", G5.x, G5.y, G5.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", H6.x, H6.y, H6.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", I7.x, I7.y, I7.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", D8.x, D8.y, D8.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", E9.x, E9.y, E9.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", F10.x, F10.y, F10.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", O11.x, O11.y, O11.z); File << (cc);

						sprintf_s(cc, "f %d %d %d %d\n", 1 + r, 2 + r, 5 + r, 3 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d %d\n", 1 + r, 3 + r, 6 + r, 4 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d %d\n", 1 + r, 4 + r, 7 + r, 2 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 5 + r, 2 + r, 8 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 5 + r, 9 + r, 3 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 6 + r, 3 + r, 9 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 6 + r, 10 + r, 4 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 7 + r, 4 + r, 10 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 7 + r, 8 + r, 2 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d %d\n", 11 + r, 9 + r, 5 + r, 8 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d %d\n", 11 + r, 10 + r, 6 + r, 9 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d %d\n", 11 + r, 8 + r, 7 + r, 10 + r); File << (cc);

						r += 11;
					}
				}
			}

			for (l = 0; h + l * z < high - h; l++) //cell2，转了60°，x坐标+0.5*x,y坐标+0.333*y,z坐标+0.5*z
			{
				for (n = 0; 0.5*sqrt(3)*a + n * y < width - 0.5*sqrt(3)*a; n++)
				{
					for (m = 0; a + m * x + n * 0.5*x < length - a; m++)
					{

						point O1 = { 0.5*a + m * x + n * 0.5*x + 0.5*x,   0.5*0.333*sqrt(3)*a + n * y + 0.333*y,   0.00 + l * z + 0.5*z };

						point D2 = { 0.00 + m * x + n * 0.5*x + 0.5*x,    a*0.33*sqrt(3) + n * y + 0.333*y,        0.25*h + l * z + 0.5*z };
						point E3 = { 0.5*a + m * x + n * 0.5*x + 0.5*x,   -0.5*0.333*sqrt(3)*a + n * y + 0.333*y,  0.25*h + l * z + 0.5*z };
						point F4 = { a + m * x + n * 0.5*x + 0.5*x,       a*0.33*sqrt(3) + n * y + 0.333*y,        0.25*h + l * z + 0.5*z };

						point G5 = { 0.00 + m * x + n * 0.5*x + 0.5*x,    0.00 + n * y + 0.333*y,                  0.5*h + l * z + 0.5*z };
						point H6 = { a + m * x + n * 0.5*x + 0.5*x,       0.00 + n * y + 0.333*y,                  0.5*h + l * z + 0.5*z };
						point I7 = { 0.5*a + m * x + n * 0.5*x + 0.5*x,   0.5*sqrt(3)*a + n * y + 0.333*y,         0.5*h + l * z + 0.5*z };

						point D8 = { 0.00 + m * x + n * 0.5*x + 0.5*x,    a*0.33*sqrt(3) + n * y + 0.333*y,        0.75*h + l * z + 0.5*z };
						point E9 = { 0.5*a + m * x + n * 0.5*x + 0.5*x,   -0.5*0.333*sqrt(3)*a + n * y + 0.333*y,  0.75*h + l * z + 0.5*z };
						point F10 = { a + m * x + n * 0.5*x + 0.5*x,       a*0.33*sqrt(3) + n * y + 0.333*y,       0.75*h + l * z + 0.5*z };

						point O11 = { 0.5*a + m * x + n * 0.5*x + 0.5*x,  0.5*0.333*sqrt(3)*a + n * y + 0.333*y,   h + l * z + 0.5*z };

						sprintf_s(cc, "v %f %f %f\n", O1.x, O1.y, O1.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", D2.x, D2.y, D2.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", E3.x, E3.y, E3.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", F4.x, F4.y, F4.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", G5.x, G5.y, G5.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", H6.x, H6.y, H6.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", I7.x, I7.y, I7.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", D8.x, D8.y, D8.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", E9.x, E9.y, E9.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", F10.x, F10.y, F10.z); File << (cc);
						sprintf_s(cc, "v %f %f %f\n", O11.x, O11.y, O11.z); File << (cc);

						sprintf_s(cc, "f %d %d %d %d\n", 1 + r, 2 + r, 5 + r, 3 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d %d\n", 1 + r, 3 + r, 6 + r, 4 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d %d\n", 1 + r, 4 + r, 7 + r, 2 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 5 + r, 2 + r, 8 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 5 + r, 9 + r, 3 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 6 + r, 3 + r, 9 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 6 + r, 10 + r, 4 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 7 + r, 4 + r, 10 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d\n", 7 + r, 8 + r, 2 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d %d\n", 11 + r, 9 + r, 5 + r, 8 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d %d\n", 11 + r, 10 + r, 6 + r, 9 + r); File << (cc);
						sprintf_s(cc, "f %d %d %d %d\n", 11 + r, 8 + r, 7 + r, 10 + r); File << (cc);

						r += 11;
					}
				}
			}


			sprintf_s(cc, "v %f %f %f\n", -thickness,         -thickness - 0.5*a*0.33*sqrt(3), high + thickness); File << (cc);
			sprintf_s(cc, "v %f %f %f\n", length + thickness, -thickness - 0.5*a*0.33*sqrt(3), high + thickness); File << (cc);
			sprintf_s(cc, "v %f %f %f\n", -thickness,         width + thickness,               high + thickness); File << (cc);
			sprintf_s(cc, "v %f %f %f\n", length + thickness, width + thickness,               high + thickness); File << (cc);

			sprintf_s(cc, "v %f %f %f\n", -thickness,         width + thickness,               -thickness - 8); File << (cc);
			sprintf_s(cc, "v %f %f %f\n", length + thickness, width + thickness,               -thickness - 8); File << (cc);
			sprintf_s(cc, "v %f %f %f\n", -thickness,         -thickness - 0.5*a*0.33*sqrt(3), -thickness - 8); File << (cc);
			sprintf_s(cc, "v %f %f %f\n", length + thickness, -thickness - 0.5*a*0.33*sqrt(3), -thickness - 8); File << (cc);

			sprintf_s(cc, "f %d %d %d %d\n", 1 + r, 2 + r, 4 + r, 3 + r); File << (cc);
			sprintf_s(cc, "f %d %d %d %d\n", 7 + r, 5 + r, 6 + r, 8 + r); File << (cc);
			sprintf_s(cc, "f %d %d %d %d\n", 1 + r, 3 + r, 5 + r, 7 + r); File << (cc);
			sprintf_s(cc, "f %d %d %d %d\n", 8 + r, 6 + r, 4 + r, 2 + r); File << (cc);
			sprintf_s(cc, "f %d %d %d %d\n", 7 + r, 8 + r, 2 + r, 1 + r); File << (cc);
			sprintf_s(cc, "f %d %d %d %d\n", 3 + r, 4 + r, 6 + r, 5 + r); File << (cc);


			File.close();
		}
		
	}


	return 0;
}
