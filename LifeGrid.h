#pragma once
#include <iostream> 
#include <vector>
#include "Kernel.h"

class LifeGrid
{
private:
	LInt *_data;
	LInt *_new_data;

	LInt *_bit0;
	LInt *_bit1;

	LInt *_single;

	int _N;
	std::vector<int> _minx, _maxx;
	std::vector<int> minx, maxx;

	LInt zero;

	void Add(LInt& b1, LInt &b0, const LInt& val);
	void Add(LInt& b2, LInt& b1, LInt &b0, const LInt& val);
	bool inline IsZero(const LInt& v);

public:

	void Init(int N);
	void Set(int x, int y, const LInt& val);
	void Set(int x, int y, int l);

	LInt Get(int x, int y);
	bool Get(int x, int y, int l);

	void ClearLevel(int l);

	void RecalculateBorders();
	void Iterate();
	void Iterate(int num_iter);
	void PrintLevel(int l);
};
