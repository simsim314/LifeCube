#include <ctime>
#include <omp.h> 
#include "GolUtils.h"

int main()
{
	std::clock_t start;
	double duration;
	
	start = std::clock();
	omp_set_num_threads(8); 
	
	#pragma omp parallel
	{
		__m256i val = _mm256_setr_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 << 8, 7 << 7, 1 << 7, 0, 0, 0);

		#pragma omp for
		for (long long i = 0; i < 10000000000ULL; i++)
			val = Iterate(val);
	
		if(omp_get_thread_num() == 0)
			print(val);
	}
	
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	
	std::cout << "duration: " << duration << '\n';
	std::cout << "BCO/s " << (10000000000.0 * 14.0 * 14.0 / 1000000000.0) / duration << '\n';

	getchar();
	return 0;
}