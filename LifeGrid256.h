#pragma once
#include <iostream> 
#include <vector>

class LifeGrid256
{
private:
	__m256i *_data;
	__m256i *_new_data;

	__m256i *_bit0;
	__m256i *_bit1;

	__m256i *_single;

	int _N;
	std::vector<int> _minx, _maxx;
	std::vector<int> minx, maxx;

	__m256i zero;

	void Add(__m256i& b1, __m256i &b0, const __m256i& val);
	void Add(__m256i& b2, __m256i& b1, __m256i &b0, const __m256i& val);
	bool inline IsZero(const __m256i& v);

public:

	void Init(int N);
	void Set(int x, int y, const __m256i& val);
	void Set(int x, int y, int l);

	__m256i Get(int x, int y);
	bool Get(int x, int y, int l);
	
	void ClearLevel(int l);

	void RecalculateBorders();
	void Iterate();
	void Iterate(int num_iter);
	void PrintLevel(int l);
};
