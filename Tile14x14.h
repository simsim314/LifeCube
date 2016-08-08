#pragma once
#include <immintrin.h>

class Tile14x14
{
private: 
	__m256i val; 

	__m256i prev;
	__m256i prev2;

	__m256i *evolved_default; 

public:

	Tile14x14();
	~Tile14x14();
};

