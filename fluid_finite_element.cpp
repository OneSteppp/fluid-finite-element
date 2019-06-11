#include<iostream>
#include<fstream>
#include<string>
#include<vector>
using namespace std;

int* read_mesh(double **point, vector<int> &number_boundary, vector<int> &boundary_point1, vector<int> &boundary_point2, vector <int> &element_node1, vector <int> &element_node2, vector <int> &element_node3)//读取坐标
{
	int number_point;//节点个数
	ifstream infile;
	double temp;//将无用信息赋值给temp

	infile.open("point.txt");//打开节点文件
	infile >> number_point;
	for (int i = 0; i < 2; i++)
	{
		point[i] = new double[number_point];
	}
	//二维数组point[0][j],[1][j]分别为第j个节点的x,y坐标
	for (int j = 0; j < number_point; j++)
	{
		infile >> temp >> point[0][j] >> point[1][j] >> temp;
	}
	infile.close();

	infile.open("boundary.txt");//打开边界文件
	while (!infile.eof())
	{
		int n;
		infile >> temp >> temp >> temp >> n;
		number_boundary.push_back(n);//边界编号，其中5为进出口边界
		infile >> temp >> n;
		boundary_point1.push_back(n - 1);//x坐标
		infile >> n;
		boundary_point2.push_back(n - 1);//y坐标
	}
	infile.close();

	infile.open("field.txt");//打开内部场单元文件
	while (!infile.eof())
	{
		int n;
		infile >> temp >> temp >> temp >> temp >> temp >> n;
		element_node1.push_back(n - 1);
		infile >> n;
		element_node2.push_back(n - 1);
		infile >> n;
		element_node3.push_back(n - 1);
	}//三角形单元的三个节点
	infile.close();

	int *p = new int[3];
	p[0] = number_point;//节点个数
	p[1] = number_boundary.size();//边界个数
	p[2] = element_node1.size();//单元个数

	return p;//返回数组指针，分别为节点个数，边界个数，单元个数
}

double* matrix_subtract(double *a, double *b, int &m)
{
	for (int i = 0; i < m; i++)
	{
		b[i] = a[i] - b[i];
	}
	return b;
}//矩阵相减

double Norm_2(double *a, int &m)
{
	double n = 0;
	for (int i = 0; i < m; i++)
	{
		n = n + a[i] * a[i];
	}
	n = sqrt(n);
	return n;
}//求范数

void Gauss_Seidel(double **A, double *y, double *r, int &m, double error)
{
	double *y_old;
	y_old = new double[m];
	double N_2;
	do
	{
		for (int i = 0; i < m; i++)
		{
			y_old[i] = y[i];
		}
		for (int i = 0; i < m; i++)
		{
			double sum = 0;
			for (int j = 0; j < m; j++)
			{
				if (j == i) continue;//跳过Aii
				sum = sum + A[i][j] * y[j];
			}
			y[i] = (r[i] - sum) / A[i][i];
		}
		N_2 = Norm_2(matrix_subtract(y, y_old, m), m);
		cout << N_2<<endl;
	} while (N_2> error);
	delete[] y_old;
	y_old = NULL;
}//高斯-赛德尔迭代法

int main()
{
	double **point = new double*[2];//节点坐标
	vector <int> number_boundary;//节点编号，其中5为出入口边界
	vector <int> boundary_point1, boundary_point2;//边界节点
	vector <int> element_node1, element_node2, element_node3;//单元内节点编号
	int* info=read_mesh(point, number_boundary, boundary_point1, boundary_point2, element_node1, element_node2, element_node3);
	int **element_node = new int*[info[2]];
	
	for (int i = 0; i < info[2]; i++)
	{
		element_node[i] = new int[3];
	}
	for (int i = 0; i < info[2]; i++)
	{
		element_node[i][0] = element_node1[i];
		element_node[i][1] = element_node2[i];
		element_node[i][2] = element_node3[i];
	}//第i个单元的第1,2,3个节点
	
	double **A = new double*[info[2]];
	double **B = new double*[info[2]];
	double **C = new double*[info[2]];
	double *E = new double[info[2]];

	for (int i = 0; i < info[2]; i++)
	{
		A[i] = new double[3];
		B[i] = new double[3];
		C[i] = new double[3];
	}
	for (int i = 0; i < info[2]; i++)
	{
        for (int j = 0; j < 3; j++)
		{
			A[i][j] = point[0][element_node[i][(j + 1) % 3]] * point[1][element_node[i][(j + 2) % 3]] - point[1][element_node[i][(j + 1) % 3]] * point[0][element_node[i][(j + 2) % 3]];

			B[i][j] = point[1][element_node[i][(j + 1) % 3]] - point[1][element_node[i][(j + 2) % 3]];

			C[i][j] = point[0][element_node[i][(j + 2) % 3]] - point[0][element_node[i][(j + 1) % 3]];
		}
		E[i] = (B[i][0] * C[i][1] - B[i][1] * C[i][0]);
	}

	double **K = new double*[info[0]];
	for (int i = 0; i < info[0]; i++)
	{
		K[i] = new double[info[0]];
	}//系数矩阵

	for (int i = 0; i < info[0]; i++)
	{
		for (int j = 0; j < info[0]; j++)
		{
			K[i][j] = 0;
		}
	}//初始为零

	for (int i = 0; i < info[2]; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				K[element_node[i][j]][element_node[i][k]] = K[element_node[i][j]][element_node[i][k]] + (B[i][j
				] * B[i][k] + C[i][j] * C[i][k]) / (4 * E[i]);
			}
		}
	}

	double *Phi = new double[info[0]];
	double *R = new double[info[0]];
	for (int i = 0; i < info[0]; i++)
	{
		Phi[i] = 0;
		R[i] = 0;
	}//待求量与右边项
	
	for (int i = 0; i < info[1]; i++)
	{
		if (number_boundary[i] == 5)
		{
			K[boundary_point1[i]][boundary_point1[i]] = 1e20;
			K[boundary_point2[i]][boundary_point2[i]] = 1e20;
			R[boundary_point1[i]] = 1e20*point[0][boundary_point1[i]];
			R[boundary_point2[i]] = 1e20*point[0][boundary_point2[i]];
		}			
	}//处理边界条件

	Gauss_Seidel(K, Phi, R, info[0], 1e-6);//求解

	//输出文件
	ofstream outfile;
	outfile.open("result.plt");
	outfile << "VARIABLES=" << "\"X\"" << "\"Y\"" << "\"Phi\"" << endl;
	outfile << "ZONE N=" << info[0] << ", E=" << info[2] << ", F=FEPOINT, ET=TRIANGLE" << endl;
	for (int i = 0; i < info[0]; i++)
	{
		outfile << point[0][i] << " " << point[1][i] << " " << Phi[i] << endl;
	}
	cout << endl;

	for (int i = 0; i < info[2]; i++)
	{
		outfile << element_node[i][0]+1 << " " << element_node[i][1]+1 << " " << element_node[i][2]+1 << endl;
	}
	outfile.close();

	//释放内存
	delete[] point[0];
	delete[] point[1];
	delete[] point;
	for (int i = 0; i < info[2]; i++)
	{
		delete[] element_node[i];
		element_node[i] = NULL;
		delete[] A[i];
		A[i] = NULL;
		delete[] B[i];
		B[i] = NULL;
		delete[] C[i];
		C[i] = NULL;
	}
	delete[] element_node;
	element_node = NULL;
	delete[] A;
	A = NULL;
	delete[] B;
	B = NULL;
	delete[] C;
	C = NULL;
	delete[] E;
	E = NULL;
	for (int i = 0; i < info[0]; i++)
	{
		delete[] K[i];
		K[i] = NULL;
	}
	delete[] K;
	K = NULL;
	delete[] Phi;
	Phi = NULL;
	delete[] R;
	R = NULL;

	return 0;
}
