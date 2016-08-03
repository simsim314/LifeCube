#include "LifeGrid.h"
#include <algorithm>
#include <malloc.h>

using namespace std;

void LifeGrid::Init(int N)
{
	_N = N;

	_data = (LInt*)_aligned_malloc(N * N * AllocSize, AllocBlock);
	_bit0 = (LInt*)_aligned_malloc(N * N * AllocSize, AllocBlock);
	_bit1 = (LInt*)_aligned_malloc(N * N * AllocSize, AllocBlock);
	_new_data = (LInt*)_aligned_malloc(N * N * AllocSize, AllocBlock);

	_minx = vector<int>(N);
	_maxx = vector<int>(N);

	minx = vector<int>(N);
	maxx = vector<int>(N);

	zero = set_zero();

	for (int i = 0; i < _N * _N; i++)
	{
		_data[i] = zero;
		_bit0[i] = zero;
		_bit1[i] = zero;
		_new_data[i] = zero;
	}

	RecalculateBorders();

	_single = (LInt*)_aligned_malloc(AllocBlock * AllocSize, AllocBlock);

	for (int i = 0; i < 64; i++)
	{
		#ifdef SSE
			_single[i] = set_epi64x(1ULL << i, 0);
			_single[i + 64] = set_epi64x(0, 1ULL << i);
		#endif

		#ifdef AVX
			_single[i] = set_epi64x(1ULL << i, 0, 0, 0);
			_single[i + 64] = set_epi64x(0, 1ULL << i, 0, 0);
			_single[i + 128] = set_epi64x(0, 0, 1ULL << i, 0);
			_single[i + 192] = set_epi64x(0, 0, 0, 1ULL << i);
		#endif
	}
}

void LifeGrid::Set(int x, int y, const LInt& val)
{
	_data[_N * y + x] = val;
}

void LifeGrid::Set(int x, int y, int l)
{
	_data[_N * y + x] = OR(_data[_N * y + x], _single[l]);
}

void LifeGrid::ClearLevel(int l)
{
	for (int j = 0; j <= _N - 1; j++)
	{
		for (int i = _minx[j]; i <= _maxx[j]; i++)
		{
			int idx = _N * j + i;
			_data[idx] = ANDNOT(_single[l], _data[idx]);
		}
	}


}

LInt LifeGrid::Get(int x, int y)
{
	return _data[_N * y + x];
}

bool LifeGrid::Get(int x, int y, int l)
{
	int idx = _N * y + x;
	LInt val = AND(_single[l], _data[idx]);
	return !IsZero(val);
}


void LifeGrid::RecalculateBorders()
{
	for (int j = 0; j < _N; j++)
	{
		minx[j] = 100000;
		maxx[j] = -100000;

		for (int i = 0; i < _N; i++)
		{
			if (!IsZero(_data[j * _N + i]))
			{
				minx[j] = min(minx[j], i);
				maxx[j] = max(maxx[j], i);
			}
		}

		if (minx[j] > 0)
			minx[j]--;

		if (maxx[j] < _N - 1)
			maxx[j]++;
	}


	for (int j = 0; j < _N; j++)
	{

		int min1 = 1000000, min2 = minx[j], min3 = 1000000;
		int max1 = -1000000, max2 = maxx[j], max3 = -1000000;

		if (j > 0)
		{
			min1 = minx[j - 1];
			max1 = maxx[j - 1];
		}

		if (j < _N - 1)
		{
			min3 = minx[j + 1];
			max3 = maxx[j + 1];
		}

		_minx[j] = min(min(min1, min2), min3);
		_maxx[j] = max(max(max1, max2), max3);
	}

}

void LifeGrid::Add(LInt& b1, LInt &b0, const LInt& val)
{
	if (IsZero(val))
		return;

	b1 = OR(AND(b0, val), b1);
	b0 = XOR(b0, val);
}

void LifeGrid::Add(LInt& b2, LInt& b1, LInt &b0, const LInt& val)
{
	if (IsZero(val))
		return;

	LInt t_b2 = AND(b0, val);

	b2 = OR(AND(t_b2, b1), b2);
	b1 = XOR(t_b2, b1);
	b0 = XOR(b0, val);
}

void LifeGrid::PrintLevel(int l)
{
	for (int i = 0; i < _N; i++)
	{
		for (int j = 0; j < _N; j++)
		{
			if (Get(j, i, l))
				cout << "X";
			else
				cout << ".";
		}

		cout << "\n";
	}

	cout << "\n";
}

bool inline LifeGrid::IsZero(const LInt& v)
{
	return  is_zero(v, v) != 0;
}

void LifeGrid::Iterate(int num_iter)
{
	for (int i = 0; i < num_iter; i++)
		Iterate();
}

void LifeGrid::Iterate()
{

	for (int j = 0; j < _N; j++)
	{
		for (int i = _minx[j]; i <= _maxx[j]; i++)
		{
			LInt val_L = zero;
			LInt val_R = zero;

			int idx = _N * j + i;

			if (i > 0)
				val_L = _data[idx - 1];

			if (i < _N - 1)
				val_R = _data[idx + 1];

			_bit1[idx] = zero;
			_bit0[idx] = _data[idx];

			Add(_bit1[idx], _bit0[idx], val_L);
			Add(_bit1[idx], _bit0[idx], val_R);
		}
	}

	for (int j = 0; j < _N; j++)
	{
		minx[j] = 1000000;
		maxx[j] = -1000000;

		for (int i = _minx[j]; i <= _maxx[j]; i++)
		{
			LInt b1_U = zero, b0_U = zero;
			LInt b1_D = zero, b0_D = zero;

			int idx = _N * j + i;

			LInt b1 = _bit1[idx], b0 = _bit0[idx], v = _data[idx];

			if (j > 0)
			{
				b1_U = _bit1[idx - _N];
				b0_U = _bit0[idx - _N];
			}

			if (j < _N - 1)
			{
				b1_D = _bit1[idx + _N];
				b0_D = _bit0[idx + _N];
			}

			LInt b2_U = zero;

			Add(b2_U, b1_U, b0_U, b0_D);
			Add(b2_U, b1_U, b1_D);

			LInt sum0 = b0_U;
			LInt sum1 = b1_U;
			LInt sum2 = b2_U;

			Add(sum2, sum1, sum0, b0);
			Add(sum2, sum1, b1);

			LInt  total3 = ANDNOT(sum2, AND(sum1, sum0));

			//Sub val from center 
			LInt b0_1 = XOR(b0, v);
			LInt b1_1 = XOR(b1, ANDNOT(b0, v));

			Add(b2_U, b1_U, b0_U, b0_1);
			Add(b2_U, b1_U, b1_1);

			LInt total3_1 = ANDNOT(b2_U, AND(b1_U, b0_U));

			_new_data[idx] = OR(total3_1, total3);

			if (!IsZero(_new_data[idx]))
			{
				minx[j] = min(minx[j], i);
				maxx[j] = max(maxx[j], i);
			}
		}

		minx[j] = max(0, minx[j] - 1);
		maxx[j] = min(_N - 1, maxx[j] + 1);
	}

	for (int j = 0; j < _N; j++)
	{
		int min1 = 1000000, min2 = minx[j], min3 = 1000000;
		int max1 = -1000000, max2 = maxx[j], max3 = -1000000;

		if (j > 0)
		{
			min1 = minx[j - 1];
			max1 = maxx[j - 1];

		}

		if (j < _N - 1)
		{
			min3 = minx[j + 1];
			max3 = maxx[j + 1];
		}

		_minx[j] = max(0, _minx[j] - 1);
		_maxx[j] = min(_N - 1, _maxx[j] + 1);

		for (int i = min(minx[j], _minx[j]); i <= max(maxx[j], _maxx[j]); i++)
		{
			int idx = j * _N + i;
			_data[idx] = _new_data[idx];
			_new_data[idx] = zero;
			_bit0[idx] = zero;
			_bit1[idx] = zero;
		}

		_minx[j] = min(min(min1, min2), min3);
		_maxx[j] = max(max(max1, max2), max3);

	}
}