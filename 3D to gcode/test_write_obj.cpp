#include <iostream>
#include <fstream>
#include<string>
using namespace std;
//将int转为string类型，后可以操作。

int main()
{
	ofstream File;
	for (int i = 0; i < 5; i++)
	{
		string bb = "";
		bb = to_string(static_cast<long long>(i));
		File.open("file_" + bb + ".obj");
		File << ("hello world " + bb +"asd");
		File.close();
	}
	return 0;
}