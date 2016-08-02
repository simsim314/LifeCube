#include "LifeGrid256.h"
#include <algorithm>
#include <malloc.h>

using namespace std; 

void LifeGrid256::Init(int N)
{
	_N = N;

	_data = (__m256i*)_aligned_malloc(N * N * 32, 256);
	_bit0 = (__m256i*)_aligned_malloc(N * N * 32, 256);
	_bit1 = (__m256i*)_aligned_malloc(N * N * 32, 256);
	_new_data = (__m256i*)_aligned_malloc(N * N * 32, 256);

	_minx = vector<int>(N);
	_maxx = vector<int>(N);

	minx = vector<int>(N);
	maxx = vector<int>(N);

	zero = _mm256_setzero_si256();

	for (int i = 0; i < _N * _N; i++)
	{
		_data[i] = zero;
		_bit0[i] = zero;
		_bit1[i] = zero;
		_new_data[i] = zero;
	}

	RecalculateBorders();

	_single = (__m256i*)_aligned_malloc(256 * 32, 256);

	for (int i = 0; i < 64; i++)
	{
		_single[i] = _mm256_set_epi64x(1ULL << i, 0, 0, 0);
		_single[i + 64] = _mm256_set_epi64x(0, 1ULL << i, 0, 0);
		_single[i + 128] = _mm256_set_epi64x(0, 0, 1ULL << i, 0);
		_single[i + 192] = _mm256_set_epi64x(0, 0, 0, 1ULL << i);
	}
}

void LifeGrid256::Set(int x, int y, const __m256i& val)
{
	_data[_N * y + x] = val;
}

void LifeGrid256::Set(int x, int y, int l)
{
	_data[_N * y + x] = _mm256_or_si256(_data[_N * y + x], _single[l]);
}

void LifeGrid256::ClearLevel(int l)
{
	for (int j = 0; j <= _N - 1; j++)
	{
		for (int i = _minx[j]; i <= _maxx[j]; i++)
		{
			int idx = _N * j + i;
			_data[idx] = _mm256_andnot_si256(_single[l], _data[idx]);
		}
	}

	
}

__m256i LifeGrid256::Get(int x, int y)
{
	return _data[_N * y + x];
}

bool LifeGrid256::Get(int x, int y, int l)
{
	int idx = _N * y + x;
	__m256i val = _mm256_and_si256(_single[l], _data[idx]);
	return !IsZero(val);
}


void LifeGrid256::RecalculateBorders()
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

void LifeGrid256::Add(__m256i& b1, __m256i &b0, const __m256i& val)
{
	if (IsZero(val))
		return;

	b1 = _mm256_or_si256(_mm256_and_si256(b0, val), b1);
	b0 = _mm256_xor_si256(b0, val);
}

void LifeGrid256::Add(__m256i& b2, __m256i& b1, __m256i &b0, const __m256i& val)
{
	if (IsZero(val))
		return;

	__m256i t_b0 = _mm256_xor_si256(b0, val);
	__m256i t_b2 = _mm256_and_si256(b0, val);
	__m256i t_b1 = _mm256_xor_si256(t_b2, b1);

	b2 = _mm256_or_si256(_mm256_and_si256(t_b2, b1), b2);
	b1 = t_b1;
	b0 = t_b0;
}

void LifeGrid256::PrintLevel(int l)
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

bool inline LifeGrid256::IsZero(const __m256i& v)
{
	return  _mm256_testz_si256(v, v) != 0;
}

void LifeGrid256::Iterate(int num_iter)
{
	for (int i = 0; i < num_iter; i++)
		Iterate();
}

void LifeGrid256::Iterate()
{

	for (int j = 0; j < _N; j++)
	{
		for (int i = _minx[j]; i <= _maxx[j]; i++)
		{
			__m256i val_L = zero;
			__m256i val_R = zero;

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
			__m256i b1_U = zero, b0_U = zero;
			__m256i b1_D = zero, b0_D = zero;

			int idx = _N * j + i;

			__m256i b1 = _bit1[idx], b0 = _bit0[idx], v = _data[idx];

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

			__m256i b2_U = zero;

			Add(b2_U, b1_U, b0_U, b0_D);
			Add(b2_U, b1_U, b1_D);

			__m256i sum0 = b0_U;
			__m256i sum1 = b1_U;
			__m256i sum2 = b2_U;

			Add(sum2, sum1, sum0, b0);
			Add(sum2, sum1, b1);

			__m256i  total3 = _mm256_and_si256(sum1, sum0);
			total3 = _mm256_andnot_si256(sum2, total3);

			//Sub val from center 
			__m256i b0_1 = _mm256_xor_si256(b0, v);
			__m256i b1_1 = _mm256_xor_si256(b1, _mm256_andnot_si256(b0, v));

			Add(b2_U, b1_U, b0_U, b0_1);
			Add(b2_U, b1_U, b1_1);

			__m256i total3_1 = total3;
			total3 = _mm256_and_si256(b1_U, b0_U);
			total3 = _mm256_andnot_si256(b2_U, total3);

			_new_data[idx] = _mm256_or_si256(total3_1, total3);

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