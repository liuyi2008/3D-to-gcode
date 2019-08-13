本程序的功能是将3D模型切片，然后2D平面进行填充，最终根据坐标再转化成gcode文件：

    3D model - slice - fill - gcode
    
第一层文件夹中的.cpp文件都是在编写过程中某些子功能的实验，具有参考价值

真正的主程序 main.cpp 在第二层3D to gcode 文件夹中
相同路径下还有：
  1、模型文件.ojb
  2、opengmesh的一些lib文件
  3、实现切片功能的slice.cpp slice.h
  
  
