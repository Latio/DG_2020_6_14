
#include   < iostream >
#include   < vector >
#include   < math.h >
using     namespace   std;

typedef vector   < double >   s_line;  //  用来表示一行
s_line line;

typedef vector   < s_line >   s_matrix;  //  用来表示一个矩阵
s_matrix matrix;
s_matrix mat;
int   nSize;  //  矩阵维数
int   nSign;  //  标记行列式正负
//void   outprint(s_matrix  &   _mat);
void   printstep(s_matrix  &   _mat);
int   step = 0;
void   line_add(s_matrix  &   _mat, int   a, int   b, double   k = 1.0)  //  第b行乘k加到第a行
{
	int   size = _mat[0].size();

	for (int i = 0; i < size; ++i)
	{
		_mat[a][i] += _mat[b][i] * k;

	}  //  end for
}



void   work1(s_matrix  &   _mat)  //  主计算函数
{

	for (int i = 1; i < nSize; ++i)
	{

		if (fabs(_mat[i - 1][i - 1]) < 0.000001)
		{
			int   mm;
			for (mm = i; mm < nSize; ++mm)
			{
				if (fabs(_mat[mm - 1][i - 1]) > 0.000001)   break;
			}  //  end for
			line_add(_mat, i - 1, mm - 1);
		}  //  end if

		for (int j = i; j < nSize; ++j)
		{
			line_add(_mat, j, i - 1, -_mat[j][i - 1] / _mat[i - 1][i - 1]);

		}  //  end for j
		//printstep(_mat);
	}  //  end for i

}


void   work2(s_matrix  &   _mat)  //  第二部计算
{
	for (int i = nSize - 2; i >= 0; --i)
	{
		for (int j = i; j >= 0; --j)
		{
			line_add(_mat, j, i + 1, -_mat[j][i + 1] / _mat[i + 1][i + 1]);
		}
		//printstep(_mat);
	}

}


void   makeunit(s_matrix  &   _mat)  //  单位化
{

	mat.clear();

	for (int i = 0; i < nSize; ++i)
	{
		line.clear();
		for (int j = 0; j < nSize * 2; ++j)
		{
			double   tmp = _mat[i][j] / _mat[i][i];
			if (fabs(tmp) < 0.000001) tmp = 0;
			line.push_back(tmp);
		}
		mat.push_back(line);
		//  cout<<endl;
	}
	_mat = mat;
}

void   printstep(s_matrix  &   _mat)  //  显示求的过程
{
	cout << "  第   " << ++step << "  步  " << endl;
	for (int i = 0; i < nSize; ++i)
	{

		for (int j = 0; j < 2 * nSize; ++j)
		{
			if (fabs(_mat[i][j]) < 0.000001) _mat[i][j] = 0;
			cout << _mat[i][j] << "     ";
			if (j == nSize - 1)cout << "   |   ";
		}
		cout << endl;
	}
	cout << endl;

}

void   outprint(s_matrix  &   _mat, double *result)  //  输出函数
{
	for (int i = 0; i < nSize; ++i)
	{

		for (int j = nSize; j < 2 * nSize; ++j)
		{
			//cout << _mat[i][j] << "     ";
			result[i + (j - nSize)*nSize] = _mat[i][j];
		}
		//cout << endl;
	}


}

void  GetMatrixInverse(int np, double *array, double *result)
{

	nSize = np;

	double *test = array;

	for (int i = 0; i < nSize; ++i)
	{
		line.clear();
		/*cout << "  输入第  " << i + 1 << "   行:   " << endl;*/

		for (int j = 0; j < nSize; ++j)
		{
			double    tmp;
			/*	cin >> tmp;*/
			tmp = test[i + j * nSize];
			line.push_back(tmp);   //  压入一个数到某行}

		}

		for (int j = 0; j < nSize; ++j)
		{
			if (i == j) line.push_back(1.0);
			else   line.push_back(0.0);
		}


		matrix.push_back(line);   //  压入一行到矩阵
	};




	work1(matrix);
	work2(matrix);
	makeunit(matrix);
	//cout << endl << "  ########################  " << endl
		//<< "  求逆结果：  " << endl;
	outprint(matrix, result);
	//cout << "  ########################  " << endl;



	//return  0;
}

//int main() {
//
//	double array[9] = { 1,3,1,2,2,4,3,1,2 };
//	double *result = new double[9];
//
//	GetMatrixInverse(3, array, result);
//
//
//	for (size_t i = 0; i < 9; i++)
//	{
//		cout << result[i] << endl;
//	}
//	delete[]result;
//
//	return 0;
//}