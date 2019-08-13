#include"slice.h"


vector<point>coord;
point p;

const float ERR = 0.001;

//�ļ���ȡ�йص�
MyMesh mesh;

extern const string file_1 = "cow_file.obj";
extern const string file_2 = "cactus.ply";
extern const string file_3 = "Armadillo.off";
int currentfile = 1;

GLuint showFaceList, showWireList, showCutList;
int showstate = 1;
bool showFace = true;
bool showWire = false;
bool showCut = true;
bool showFlatlines = false;

//���ڰ�Χ�������
double Bmax_x, Bmax_y, Bmax_z, Bmin_x, Bmin_y, Bmin_z, px, py, pz;
const int INF_formistake = 0x3f3f3f3f;
int mid_n = 0; int mid_m = 0;
int countt = 0;

MyMesh::Point plist[3000][1000];
int numofcut[3000];

// ��ȡ�ļ��ĺ���
void readfile(string file) {
	// ���󶥵㷨�� vertex normals
	mesh.request_vertex_normals();
	//��������ڶ��㷨�ߣ��򱨴� 
	if (!mesh.has_vertex_normals())
	{
		cout << "���󣺱�׼�������� �����ߡ�������" << endl;
		return;
	}
	// ����ж��㷢�����ȡ�ļ� 
	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, file, opt))
	{
		cout << "�޷���ȡ�ļ�:" << file << endl;
		return;
	}
	else cout << "�ɹ���ȡ�ļ�:" << file << endl;
	cout << endl; // Ϊ��ui��ʾ�ÿ�һЩ
				  //��������ڶ��㷨�ߣ�������
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		// ͨ���淨�߼��㶥�㷨��
		mesh.request_face_normals();
		// mesh��������㷨��
		mesh.update_normals();
		// �ͷ��淨��
		mesh.release_face_normals();
	}
}

//ͨ���򵥵ı�������ģ�Ͱ�Χ�У����ȡ�����޹أ�
void BoundingBox() {
	MyMesh::Point pt;
	int st = 0;
	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
		MyMesh::VertexHandle vh_i = *it;
		pt = mesh.point(vh_i);
		px = pt.data()[0];
		py = pt.data()[1];
		pz = pt.data()[2];
		if (st == 0) {
			Bmax_x = Bmin_x = px;
			Bmax_y = Bmin_y = py;
			Bmax_z = Bmin_z = pz;
			st++;
		}
		else {
			if (px > Bmax_x)Bmax_x = px; else if (px < Bmin_x)Bmin_x = px;
			if (py > Bmax_y)Bmax_y = py; else if (py < Bmin_y)Bmin_y = py;
			if (pz > Bmax_z)Bmax_z = pz; else if (pz < Bmin_z)Bmin_z = pz;
		}
	}
	//	Bmax_x += 0.2; Bmax_y += 0.2; Bmax_z += 0.2;
	//	Bmin_x -= 0.2; Bmin_y -= 0.2; Bmin_z -= 0.2;
	//printf("%f %f %f %f %f %f\n", Bmax_x, Bmax_y, Bmax_z, Bmin_x, Bmin_y, Bmin_z);
}

//���溯��
bool IntersectPlane(MyMesh::Point pt, MyMesh::Point pnorm, MyMesh::Point *pilist, int &pnum) {
	//�������� pt��pnorm��*pilist��pnum[]  ���庯��ԭ�� �� �����㷨.docx
	int starte, ne, ne1, nf, num;
	MyMesh::Point vt1, vt2;
	//MyMesh::Face f1;
	MyMesh::HalfedgeHandle nhe;
	MyMesh::FaceHandle nhf;
	float d1, d2, sd1, sd2;
	bool *flag, suc;
	float dist, mind = 1.0e+8;

	pnum = 0;
	sd1 = sd2 = -10000;
	int esize = mesh.n_halfedges();
	flag = new bool[esize];

	suc = false;


	for (MyMesh::HalfedgeIter it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) {

		MyMesh::HalfedgeHandle hh = *it;
		int id = hh.idx();
		flag[id] = false;

		auto fromVertex = mesh.from_vertex_handle(hh);
		auto toVertex = mesh.to_vertex_handle(hh);
		vt1 = mesh.point(fromVertex);
		vt2 = mesh.point(toVertex);
		//printf("$ %.3f %.3f $\n", vt1.data()[0],vt2.data()[0]);
		d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
		d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);

		if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0) || d1 > 0 && d2 == 0 || d1 == 0 && d2 > 0))
		{
			flag[id] = true;

			vt1.data()[0] = vt1.data()[0] - pt.data()[0];
			vt1.data()[1] = vt1.data()[1] - pt.data()[1];
			vt1.data()[2] = vt1.data()[2] - pt.data()[2];       // point date minus point date 
			dist = vt1.data()[0] * vt1.data()[0] + vt1.data()[1] * vt1.data()[1] + vt1.data()[2] * vt1.data()[2];
			if (dist < mind) {
				nhe = hh;
				mind = dist;
				ne = id;                    //  printf("ne:  %d  \n", ne);
				suc = true;
			}
		}
	}

	if (!suc) {
		delete[]flag;
		return false;
	}

	starte = ne;
	num = 0;

	suc = false;

	nhf = mesh.face_handle(nhe);

	while (!suc) {
		//printf("%%%%");	

		auto fromVertex = mesh.from_vertex_handle(nhe);
		auto toVertex = mesh.to_vertex_handle(nhe);

		vt1 = mesh.point(fromVertex);
		vt2 = mesh.point(toVertex);

		d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
		d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);
		//printf("$$$%lf %lf \n", d1, d2);
		if ((sd1 == d1) && (sd2 == d2)) {
			flag[ne] = false;
		}
		sd1 = d1; sd2 = d2;
		pilist[num].data()[0] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0];
		pilist[num].data()[1] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1];
		pilist[num].data()[2] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[2] - vt1.data()[2]) + vt1.data()[2];
		//printf("$$$%lf %lf %lf %lf %lf %lf\n", vt1.data()[0], vt1.data()[1], vt1.data()[2], vt2.data()[0], vt2.data()[1], vt2.data()[2]);
		num++;



		int nn = 0;
		do {
			for (auto it = mesh.fh_begin(nhf); it != mesh.fh_end(nhf); ++it) {

				MyMesh::HalfedgeHandle halfnow = *it;
				nn++;
				if (nn > 3) {
					starte = ne;
					flag[ne] = false;
					break;
				}
				const int ne1 = halfnow.idx();

				if (flag[ne1] == false || ne == ne1) continue;
				nn = 0;
				MyMesh::VertexHandle fromV = mesh.from_vertex_handle(halfnow);
				MyMesh::VertexHandle toV = mesh.to_vertex_handle(halfnow);

				vt1 = mesh.point(fromV);
				vt2 = mesh.point(toV);

				d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
					+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
				d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
					+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);

				if (((d1 >= 0 && d2 <= 0) || (d1 <= 0 && d2 >= 0)) && fabs(d1 - d2) > ERR) //ERR��ֵδ֪ 
				{

					MyMesh::HalfedgeHandle halfnext = mesh.opposite_halfedge_handle(halfnow);

					nhf = mesh.face_handle(halfnext);


					int ne2 = halfnext.idx();

					flag[ne1] = flag[ne2] = false;


					if (nhf.idx() == -1)
					{
						starte = ne;
						flag[ne] = false;
						break;
					}
					ne = ne2;
					// compute the intersecting point		
					pilist[num].data()[0] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0];
					pilist[num].data()[1] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1];
					pilist[num].data()[2] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[2] - vt1.data()[2]) + vt1.data()[2];
					//printf("##%lf %lf %lf\n", pilist[num].data()[0], pilist[num].data()[1], pilist[num].data()[2]);
					num++;

					break;
				}
			}
		} while (ne != starte);

		suc = true;

		for (auto it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) {

			MyMesh::HalfedgeHandle hh = *it;

			int id = hh.idx();

			if (flag[id] == true) {

				ne = id;

				starte = ne;
				nhe = hh;
				nhf = mesh.face_handle(nhe);
				if (nhf.idx() == -1)
				{
					flag[ne] = false;

					continue;
				}
				pilist[num].data()[0] = -10000;
				pilist[num].data()[1] = -10000;
				pilist[num].data()[2] = -10000;
				num++;


				suc = false;
				break;
			}
		}


	};

	pnum = num;
	//printf("%d\n", pnum);

	delete[]flag;

	return true;
}

//bool IntersectPlane(MyMesh::Point pt, MyMesh::Point pnorm, MyMesh::Point *pilist, int &pnum);
void findIntersect() {
	countt = 0;
	/*printf("\n\n\n--------------------------------------------------------------------------------------------\n");
	printf("          -----------------------------------------------------------------------\n\n");
	printf("             **************         ***          ***          *****************         \n");
	printf("             **************         ***          ***          *****************         \n");
	printf("             ***                    ***          ***                *****               \n");
	printf("             ***                    ***          ***                *****               \n");
	printf("             ***                    ***          ***                *****               \n");
	printf("             ***                    ***          ***                *****               \n");
	printf("             ***                    ***          ***                *****               \n");
	printf("             ***                    ***          ***                *****               \n");
	printf("             ***                    ***          ***                *****               \n");
	printf("             **************         ****************                *****               \n");
	printf("             **************         ****************                *****               \n");
	printf("          -----------------------------------------------------------------------\n");
	printf("\n--------------------------------------------------------------------------------------------\n\n");
	cout << "Note 1: From the left side of the model coordinate system (left-handed coordinate system)" << endl;
	cout << "to the right side of the section, the case is from the head to the tail of the cow model." << endl;
	printf("\n");
	cout << "Note 2: 'EX' represents an invalid data point,Often occurs when a section contains multiple section rings." << endl;
	cout << "'XX%' represents the percentage of the lengthof the cow's body," << endl;
	cout << "'part' represents which part of the cow may be taken from this section." << endl;
	printf("\n--------------------------------------------------------------------------------------------\n");
	printf("--------------------------------------------------------------------------------------------\n\n\n");*/
	printf("��Ƭ��\n");
	for (double i = Bmin_x; i < Bmax_x; ) {
		MyMesh::Point pt;
		MyMesh::Point pilist[1000];
		//MyMesh::Normal v0(0, 1, 1), v1(0, 1, -1);
		//MyMesh::Normal vf = v0 % v1;
		MyMesh::Normal vf(1, 0, 0);
		int pnum = 0;
		float Xport = 0;
		pt.data()[0] = i; pt.data()[1] = pt.data()[2] = 0;
		if (IntersectPlane(pt, vf, pilist, pnum)) {
			numofcut[countt] = pnum;
			if (Bmin_x * Bmax_x < 0) {

				if (i < 0) {
					Xport = 100 * (abs(i - Bmin_x) / (abs(Bmax_x) + abs(Bmin_x)));
				}
				else {
					Xport = 100 * ((i + abs(Bmin_x)) / (abs(Bmax_x) + abs(Bmin_x)));
				}
			}
			else {
				Xport = 100 * (abs(i - Bmin_x) / abs(abs(Bmax_x) - abs(Bmin_x)));
			}
			//cout << "######"<< Xport << endl;
			//printf("-------------------------------------------------\n");
			if (Xport < 18) {
				//printf("��%d�����滷\t    %.4f%%    ͷ��\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport);
			}
			else if (Xport < 25) {
				//printf("��%d�����滷\t    %.4f%%    ����\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport);
			}
			else if (Xport < 35) {
				//printf("��%d�����滷\t    %.4f%%    ǰ��\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport);
			}
			else if (Xport < 77) {
				//printf("��%d�����滷\t    %.4f%%    ����\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport);
			}
			else if (Xport < 96) {
				//printf("��%d�����滷\t    %.4f%%    �β�\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport);
			}
			else {
				//printf("��%d�����滷\t    %.4f%%    β��\n", countt + 1, Xport);
				printf("%.4f%%\r", Xport);
			}
			int j; int numdiff = 0; MyMesh::Point diff = pilist[0];
			for (j = 0; j < pnum; j++) {
				plist[countt][j] = pilist[j];
				if (j == pnum - 1 || ((pilist[j + 1].data()[0] == -10000) && (pilist[j + 3].data()[0] != -10000))) {

					if ((pilist[j].data()[2] - diff.data()[2]) < 0.001) {
						numdiff++;
						continue;
					}//cout << diff.data()[2] << pilist[j].data()[2] << endl;
				}
				if (pilist[j].data()[0] == -10000) {
					if (j < pnum - 1)diff = pilist[j + 1];
					//printf("-- %d th:(   EX  , EX  , EX   ) \t\t", j + 1 - numdiff);
				}
				else
				{
					//printf("-- %d th:( %.4f, %.4f, %.4f )   \t", j + 1 - numdiff, pilist[j].data()[0], pilist[j].data()[1], pilist[j].data()[2]);
					p = { pilist[j].data()[1] ,pilist[j].data()[2] };
					coord.push_back(p);
				}
				if ((j + numdiff) % 2) {
					//printf("\n");
				}
			}

			/*
			int j;
			for (j = 0; j < pnum; j++) {
				plist[countt][j] = pilist[j];
				if (pilist[j].data()[0] == -10000) {
					printf("-- %d th:(   EX  , EX  , EX   ) \t\t",j+1);
				}else printf("-- %d th:( %.4f, %.4f, %.4f )   \t", j+1, pilist[j].data()[0], pilist[j].data()[1], pilist[j].data()[2]);
				if (j % 2) {
					printf("\n");
				}
			}
			*/
			if (j % 2) {
				//printf("\n");
			}
			//printf("-------------------------------------------------\n");
			model.push_back(coord);
			coord.clear();
			countt++;
		}
		/*
		//ѡ���ԵĶԹؼ���λ��߽����ܶ�
		if ((i > Bmin_x + abs(Bmax_x - Bmin_x)*0.75 && i < Bmin_x + abs(Bmax_x - Bmin_x)*0.87) ||
			(i<Bmin_x + abs(Bmax_x - Bmin_x)*0.45 && i>Bmin_x + abs(Bmax_x - Bmin_x)*0.33) ||
			(i > Bmin_x + abs(Bmax_x - Bmin_x)*0.92)) {
			i += 20;
		}
		else {
			i += 5;
		}
		*/
		i += 5;
	}
}
