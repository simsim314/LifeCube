#include <iostream>
#include <immintrin.h>

using namespace std;

void print(const unsigned short& val)
{
	unsigned short t = val;
	bool vals[16];

	for (int i = 0; i < 16; i++)
	{
		vals[15 - i] = t % 2 == 1;
		t = t / 2;
	}

	for (int i = 0; i < 16; i++)
	{
		if (vals[i])
			cout << "O";
		else
			cout << ".";
	}

	cout << "\n";
}

void print(const __m256i& val)
{
	unsigned short* f = (unsigned short*)&val;

	for (int i = 0; i < 16; i++)
		print(f[i]);

	cout << "\n\n";
}

void Add(__m256i& b1, __m256i &b0, const __m256i& val)
{
	b1 = _mm256_or_si256(_mm256_and_si256(b0, val), b1);
	b0 = _mm256_xor_si256(b0, val);
}

void Add_Init(__m256i& b1, __m256i &b0, const __m256i& val)
{
	b1 = _mm256_and_si256(b0, val);
	b0 = _mm256_xor_si256(b0, val);
}

void Add(__m256i& b2, __m256i& b1, __m256i &b0, const __m256i& val)
{
	__m256i t_b2 = _mm256_and_si256(b0, val);

	b2 = _mm256_or_si256(_mm256_and_si256(t_b2, b1), b2);
	b1 = _mm256_xor_si256(t_b2, b1);
	b0 = _mm256_xor_si256(b0, val);
}

void Add_Init(__m256i& b2, __m256i& b1, __m256i &b0, const __m256i& val)
{
	__m256i t_b2 = _mm256_and_si256(b0, val);

	b2 = _mm256_and_si256(t_b2, b1);
	b1 = _mm256_xor_si256(t_b2, b1);
	b0 = _mm256_xor_si256(b0, val);
}

__m256i Iterate(const __m256i& val)
{
	__m256i sum_bit_0 = _mm256_slli_epi16(val, 1);
	__m256i sum_bit_1;
	__m256i sum_at_least_4;

	__m256i next_cell = _mm256_srai_epi16(val, 1);
	Add_Init(sum_bit_1, sum_bit_0, next_cell);

	const int shuffle_up = _MM_SHUFFLE(2, 0, 0, 1);
	const int shuffle_lo = _MM_SHUFFLE(0, 0, 2, 0);

	__m256i upper_word = _mm256_alignr_epi8(_mm256_permute2x128_si256(val, val, shuffle_up), val, 2);
	__m256i lower_word = _mm256_alignr_epi8(val, _mm256_permute2x128_si256(val, val, shuffle_lo), 14);

	Add(sum_bit_1, sum_bit_0, upper_word);

	next_cell = _mm256_slli_epi16(upper_word, 1);
	Add_Init(sum_at_least_4, sum_bit_1, sum_bit_0, next_cell);

	next_cell = _mm256_srai_epi16(upper_word, 1);
	Add(sum_at_least_4, sum_bit_1, sum_bit_0, next_cell);

	__m256i l_sum_bit_0 = _mm256_srai_epi16(lower_word, 1);
	__m256i l_sum_bit_1;

	Add_Init(l_sum_bit_1, l_sum_bit_0, lower_word);

	next_cell = _mm256_slli_epi16(lower_word, 1);
	Add(l_sum_bit_1, l_sum_bit_0, next_cell);

	Add(sum_at_least_4, sum_bit_1, sum_bit_0, l_sum_bit_0);
	Add(sum_at_least_4, sum_bit_1, l_sum_bit_1);

	sum_bit_0 = _mm256_or_si256(sum_bit_0, val);
	sum_bit_1 = _mm256_and_si256(sum_bit_1, sum_bit_0);
	sum_bit_1 = _mm256_andnot_si256(sum_at_least_4, sum_bit_1);

	return sum_bit_1;
}

void CopyLR(__m256i& left, __m256i& right)
{
	const __m256i maskLeft = _mm256_setr_epi16(0, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 0);
	const __m256i maskRight= _mm256_setr_epi16(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0);
	
	left = _mm256_or_si256(left, _mm256_and_si256(maskRight, _mm256_srai_epi16(right, 14)));
	right = _mm256_or_si256(right, _mm256_and_si256(maskLeft, _mm256_slli_epi16(left, 14)));
}

void CopyULBR(__m256i& upper_left, __m256i& bottom_right)
{
	const int shuff_up = _MM_SHUFFLE(0, 0, 2, 0);
	const int shuff_bo = _MM_SHUFFLE(2, 0, 0, 1);

	const __m256i blend_UL = _mm256_setr_epi32(0, 0, 0, 0, 0, 0, 0, 1);
	const __m256i blend_BR = _mm256_setr_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32768, 0);
	
	__m256i ul = _mm256_slli_si256(_mm256_permute2x128_si256(upper_left, upper_left, shuff_up), 13);
	ul = _mm256_slli_epi16(ul, 7);

	bottom_right = _mm256_blendv_epi8(bottom_right, ul, blend_BR);

	__m256i br = _mm256_or_si256(bottom_right, _mm256_srli_si256(_mm256_permute2x128_si256(bottom_right, bottom_right, shuff_bo), 13));
	br = _mm256_srai_epi16(br, 7);

	upper_left = _mm256_blendv_epi8(upper_left, br, blend_UL);
}


void CopyURBL(__m256i& upper_right, __m256i& bottom_left)
{
	const int shuff_up = _MM_SHUFFLE(0, 0, 2, 0);
	const int shuff_bo = _MM_SHUFFLE(2, 0, 0, 1);

	const __m256i blend_UL = _mm256_setr_epi32(0, 0, 0, 0, 0, 0, 0, 2);
	const __m256i blend_BR = _mm256_setr_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16384, 0);

	__m256i ul = _mm256_slli_si256(_mm256_permute2x128_si256(upper_right, upper_right, shuff_up), 13);
	ul = _mm256_srai_epi16(ul, 7);

	bottom_left = _mm256_blendv_epi8(bottom_left, ul, blend_BR);

	__m256i br = _mm256_or_si256(bottom_left, _mm256_srli_si256(_mm256_permute2x128_si256(bottom_left, bottom_left, shuff_bo), 13));
	br = _mm256_slli_epi16(br, 7);

	upper_right = _mm256_blendv_epi8(upper_right, br, blend_UL);
}
