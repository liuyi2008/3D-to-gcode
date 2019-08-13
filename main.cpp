#include"slice.h"

extern const string file_1;
vector<vector<point> >model;

#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
typedef unsigned char boolean;

vector<point>polypoint = { {250,50},{550,150},{550,400},{250,250},{100,350},{100,100},{120,30} }; //���������꼯
vector<point>tripoint;  //�仯�� �߶����� �ɶ�
vector<point>new_tripoint;  //�±仯�� �߶����� �ɶ�
vector<point>interpoint; //��������

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

int main(int argc, char** argv)
{

	readfile(file_1);
	BoundingBox();
	findIntersect();

	for (int i = 0; i < model.size(); i++)
	{
		for (int j = 0; j < model[i].size(); j++)
		{
			printf("x = %d, y = %f,z = %f",i, model[i][j].x, model[i][j].y);
		}

	}

	return 0;
}