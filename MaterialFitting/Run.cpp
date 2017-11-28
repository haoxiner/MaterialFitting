#include "CalculateGridWeight.h"
#include <string>
#include <Windows.h>
#undef max
std::vector<hx::Float> GenAlpha()
{
	hx::Float alphaMin = 0.007;
	hx::Float alphaMax = 0.4;
	auto logAMin = std::log(alphaMin);
	auto logAMax = std::log(alphaMax);
	auto diff = logAMax - logAMin;
	std::vector<hx::Float> ret;
	for (int i = 0; i < 40; i++) {
		auto a = std::exp(diff * i / 39 + logAMin);
		ret.push_back(a);
	}
	return ret;
}

void TestGenAlpha()
{
	auto as = GenAlpha();
	int cnt = 1;
	for (auto v : as) {
		std::cerr << cnt << ": " << v << std::endl;
		++cnt;
	}
}

void GenGridPdf()
{
	std::ifstream file(R"(F:/gridSample.bin)", std::ios::binary);
	std::ofstream out(R"(F:/gridWeight.bin)", std::ios::binary);
	std::vector<unsigned int> src(NUM_OF_GRID);
	file.read((char*)src.data(), src.size() * sizeof(unsigned int));
	file.close();
	std::vector<hx::Float> dest(NUM_OF_GRID);
	for (int i = 0; i < src.size(); i++) {
		dest[i] = src[i];
		dest[i] /= NUM_OF_GRID;
		dest[i] /= SAMPLES_PER_GRID;
	}
	out.write((char*)dest.data(), dest.size() * sizeof(hx::Float));
	out.close();
}

void GenGridCenter()
{
	std::vector<hx::Float3> gridCenter(NUM_OF_GRID);
	for (int y = 0; y < WIDTH; y++) {
		for (int x = 0; x < WIDTH; x++) {
			hx::Float px = (x + 0.5) / WIDTH * SQRT2;
			hx::Float py = (y + 0.5) / WIDTH * SQRT2;
			gridCenter[x + y * WIDTH] = hx::MapPlaneToUpperSphere({ SIN_45_DEGREE * (px + py) - 1.0, SIN_45_DEGREE * (py - px) });
		}
	}
	std::ofstream file(R"(F:/gridCenter.bin)", std::ios::binary);
	file.write((char*)gridCenter.data(), gridCenter.size() * sizeof(gridCenter[0]));
	file.close();
}

void GenVMF_WARD()
{
	// grid weight
	std::ifstream gridWeightFile(R"(F:/gridWeight.bin)", std::ios::binary);
	std::vector<hx::Float> gridWeight(NUM_OF_GRID);
	gridWeightFile.read((char*)gridWeight.data(), gridWeight.size() * sizeof(gridWeight[0]));
	gridWeightFile.close();

	// grid center
	std::vector<hx::Float3> gridCenter(NUM_OF_GRID);
	std::ifstream gridCenterFile(R"(F:/gridCenter.bin)", std::ios::binary);
	gridCenterFile.read((char*)gridCenter.data(), gridCenter.size() * sizeof(gridCenter[0]));
	gridCenterFile.close();

	
	std::ofstream wo(R"(F:/wo.bin)", std::ios::binary);
	const int NUM_OF_DEGREE_SAMPLES = 1024;
	for (int degreeSampleIdx = 0; degreeSampleIdx < NUM_OF_DEGREE_SAMPLES; degreeSampleIdx++) {
		hx::Float degree = degreeSampleIdx / (double)NUM_OF_DEGREE_SAMPLES * 90.0;
		auto radian = hx::DegreeToRadian(static_cast<hx::Float>(degree));
		wo.write((char*)&radian, sizeof(radian));
	}
	wo.close();
	// grid ward pdf
	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point beginning = myclock::now();
	// obtain a seed from the timer
	myclock::duration dt = myclock::now() - beginning;
	unsigned seed = dt.count();
	std::vector<std::default_random_engine> generators(NUM_OF_THREAD);
	for (int i = 0; i < NUM_OF_THREAD; i++) {
		generators[i].seed(seed + i);
	}

	std::vector<hx::Float> alphas = GenAlpha();
	omp_set_num_threads(NUM_OF_THREAD);
#pragma omp parallel for
	for (int alphaIdx = 0; alphaIdx < alphas.size(); alphaIdx++) {
		std::uniform_real_distribution<hx::Float> uniformDistributionInclusive(0.0, std::nextafter(1.0, std::numeric_limits<hx::Float>::max()));
		std::uniform_real_distribution<hx::Float> uniformDistributionExclusive(0.0, 1.0);
		auto threadID = omp_get_thread_num();
		auto& generator = generators[threadID];
		std::vector<hx::Float> gridWardPdf(NUM_OF_GRID);
		auto alpha = alphas[alphaIdx];
		std::ofstream alphaFile("F:/alpha" + std::to_string(alphaIdx + 1) + ".bin", std::ios::binary);
		for (int degreeSampleIdx = 0; degreeSampleIdx < NUM_OF_DEGREE_SAMPLES; degreeSampleIdx++) {
			hx::Float degree = degreeSampleIdx / (double) NUM_OF_DEGREE_SAMPLES * 90.0;
			hx::Float cosTheta = std::cos(hx::DegreeToRadian(static_cast<hx::Float>(degree)));
			hx::Float3 wo = { std::sqrt(1.0 - cosTheta * cosTheta), 0.0, cosTheta };

			std::vector<hx::Float> pdfX(WIDTH);
			for (int y = 0; y < WIDTH; y++) {
				for (int x = 0; x < WIDTH; x++) {
					int idx = x + y * WIDTH;
					//hx::Float px = (x + 0.5) / WIDTH * SQRT2;
					//hx::Float py = (y + 0.5) / WIDTH * SQRT2;
					//auto gridCenter = hx::MapPlaneToUpperSphere({ SIN_45_DEGREE * (px + py) - 1.0, SIN_45_DEGREE * (py - px) });

					gridWardPdf[idx] = hx::Ward(alpha, gridCenter[idx], wo) * gridWeight[idx];
					pdfX[x] += gridWardPdf[idx];
				}
			}

			// sum of pdf must be 1.0
			auto sum = std::accumulate(gridWardPdf.begin(), gridWardPdf.end(), 0.0);
			for (auto&& v : gridWardPdf) {
				v /= sum;
			}
			for (auto&& v : pdfX) {
				v /= sum;
			}
			// cdf of x
			std::vector<hx::Float> cdfX(WIDTH);
			cdfX[0] = pdfX[0];
			for (int i = 1; i < WIDTH; i++) {
				cdfX[i] = cdfX[i - 1] + pdfX[i];
			}
			// inverse cdf
			std::vector<hx::Float> cdfYAfterX(WIDTH);

			// sample vmf
			hx::Float3 vmfVectorSum;
			const int NUM_OF_SAMPLES = 4096 * 4;
			for (int sampleIdx = 0; sampleIdx < NUM_OF_SAMPLES; sampleIdx++) {
				auto lowp = std::lower_bound(cdfX.begin(), cdfX.end(), uniformDistributionInclusive(generator));
				int x = lowp - cdfX.begin();
				if (pdfX[x] == 0.0) {
					continue;
				}
				// cdf of y after x
				int idx = x + 0 * WIDTH;
				cdfYAfterX[0] = gridWardPdf[idx] / pdfX[x];
				for (int y = 1; y < WIDTH; y++) {
					idx = x + y * WIDTH;
					cdfYAfterX[y] = cdfYAfterX[y - 1] + gridWardPdf[idx] / pdfX[x];
				}
				lowp = std::lower_bound(cdfYAfterX.begin(), cdfYAfterX.end(), uniformDistributionInclusive(generator));
				int y = lowp - cdfYAfterX.begin();

				// sample vector in grid
				hx::Float px = ((hx::Float)x + uniformDistributionExclusive(generator)) / WIDTH * SQRT2;
				hx::Float py = ((hx::Float)y + uniformDistributionExclusive(generator)) / WIDTH * SQRT2;
				vmfVectorSum += hx::MapPlaneToUpperSphere({ SIN_45_DEGREE * (px + py) - 1.0, SIN_45_DEGREE * (py - px) });
			}
			hx::Float R = vmfVectorSum.Length() / NUM_OF_SAMPLES;
			auto R2 = R * R;
			auto k = R * (3 - R2) / (1 - R2);
			alphaFile.write((char*)&k, sizeof(k));
		}
	}
}

void GenVMF_GGX()
{
	// grid weight
	std::ifstream gridWeightFile(R"(F:/gridWeight.bin)", std::ios::binary);
	std::vector<hx::Float> gridWeight(NUM_OF_GRID);
	gridWeightFile.read((char*)gridWeight.data(), gridWeight.size() * sizeof(gridWeight[0]));
	gridWeightFile.close();

	// grid center
	std::vector<hx::Float3> gridCenter(NUM_OF_GRID);
	std::ifstream gridCenterFile(R"(F:/gridCenter.bin)", std::ios::binary);
	gridCenterFile.read((char*)gridCenter.data(), gridCenter.size() * sizeof(gridCenter[0]));
	gridCenterFile.close();


	std::ofstream wo(R"(F:/wo.bin)", std::ios::binary);
	const int NUM_OF_DEGREE_SAMPLES = 1024;
	for (int degreeSampleIdx = 0; degreeSampleIdx < NUM_OF_DEGREE_SAMPLES; degreeSampleIdx++) {
		hx::Float degree = degreeSampleIdx / (double)NUM_OF_DEGREE_SAMPLES * 90.0;
		auto radian = hx::DegreeToRadian(static_cast<hx::Float>(degree));
		wo.write((char*)&radian, sizeof(radian));
	}
	wo.close();
	// grid ward pdf
	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point beginning = myclock::now();
	// obtain a seed from the timer
	myclock::duration dt = myclock::now() - beginning;
	unsigned seed = dt.count();
	std::vector<std::default_random_engine> generators(NUM_OF_THREAD);
	for (int i = 0; i < NUM_OF_THREAD; i++) {
		generators[i].seed(seed + i);
	}

	std::vector<hx::Float> alphas = GenAlpha();
	omp_set_num_threads(NUM_OF_THREAD);
#pragma omp parallel for
	for (int alphaIdx = 0; alphaIdx < alphas.size(); alphaIdx++) {
		std::uniform_real_distribution<hx::Float> uniformDistributionInclusive(0.0, std::nextafter(1.0, std::numeric_limits<hx::Float>::max()));
		std::uniform_real_distribution<hx::Float> uniformDistributionExclusive(0.0, 1.0);
		auto threadID = omp_get_thread_num();
		auto& generator = generators[threadID];
		std::vector<hx::Float> gridGGXPdf(NUM_OF_GRID);
		auto alpha = alphas[alphaIdx];
		std::ofstream alphaFile("F:/ggx_alpha" + std::to_string(alphaIdx + 1) + ".bin", std::ios::binary);
		for (int degreeSampleIdx = 0; degreeSampleIdx < NUM_OF_DEGREE_SAMPLES; degreeSampleIdx++) {
			hx::Float degree = degreeSampleIdx / (double)NUM_OF_DEGREE_SAMPLES * 90.0;
			hx::Float cosTheta = std::cos(hx::DegreeToRadian(static_cast<hx::Float>(degree)));

			std::vector<hx::Float> pdfX(WIDTH);
			for (int y = 0; y < WIDTH; y++) {
				for (int x = 0; x < WIDTH; x++) {
					int idx = x + y * WIDTH;
					//hx::Float px = (x + 0.5) / WIDTH * SQRT2;
					//hx::Float py = (y + 0.5) / WIDTH * SQRT2;
					//auto gridCenter = hx::MapPlaneToUpperSphere({ SIN_45_DEGREE * (px + py) - 1.0, SIN_45_DEGREE * (py - px) });

					gridGGXPdf[idx] = hx::GGX(alpha, cosTheta) * gridWeight[idx];
					pdfX[x] += gridGGXPdf[idx];
				}
			}

			// sum of pdf must be 1.0
			auto sum = std::accumulate(gridGGXPdf.begin(), gridGGXPdf.end(), 0.0);
			for (auto&& v : gridGGXPdf) {
				v /= sum;
			}
			for (auto&& v : pdfX) {
				v /= sum;
			}
			// cdf of x
			std::vector<hx::Float> cdfX(WIDTH);
			cdfX[0] = pdfX[0];
			for (int i = 1; i < WIDTH; i++) {
				cdfX[i] = cdfX[i - 1] + pdfX[i];
			}
			// inverse cdf
			std::vector<hx::Float> cdfYAfterX(WIDTH);

			// sample vmf
			hx::Float3 vmfVectorSum;
			const int NUM_OF_SAMPLES = 4096 * 4;
			for (int sampleIdx = 0; sampleIdx < NUM_OF_SAMPLES; sampleIdx++) {
				auto lowp = std::lower_bound(cdfX.begin(), cdfX.end(), uniformDistributionInclusive(generator));
				int x = lowp - cdfX.begin();
				if (pdfX[x] == 0.0) {
					continue;
				}
				// cdf of y after x
				int idx = x + 0 * WIDTH;
				cdfYAfterX[0] = gridGGXPdf[idx] / pdfX[x];
				for (int y = 1; y < WIDTH; y++) {
					idx = x + y * WIDTH;
					cdfYAfterX[y] = cdfYAfterX[y - 1] + gridGGXPdf[idx] / pdfX[x];
				}
				lowp = std::lower_bound(cdfYAfterX.begin(), cdfYAfterX.end(), uniformDistributionInclusive(generator));
				int y = lowp - cdfYAfterX.begin();

				// sample vector in grid
				hx::Float px = ((hx::Float)x + uniformDistributionExclusive(generator)) / WIDTH * SQRT2;
				hx::Float py = ((hx::Float)y + uniformDistributionExclusive(generator)) / WIDTH * SQRT2;
				vmfVectorSum += hx::MapPlaneToUpperSphere({ SIN_45_DEGREE * (px + py) - 1.0, SIN_45_DEGREE * (py - px) });
			}
			hx::Float R = vmfVectorSum.Length() / NUM_OF_SAMPLES;
			auto R2 = R * R;
			auto k = R * (3 - R2) / (1 - R2);
			alphaFile.write((char*)&k, sizeof(k));
		}
	}
}

int main()
{
	LARGE_INTEGER t1, t2, t3, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
	std::cerr << "Running..." << std::endl;
	//CalculateGridWeight();
	//GenGridPdf();
	//GenGridCenter();
	GenVMF_WARD();
	QueryPerformanceCounter(&t2);
	std::cerr << (t2.QuadPart - t1.QuadPart) / tc.QuadPart << "s" << std::endl;

	std::cerr << "Gen GGX" << std::endl;
	GenVMF_GGX();
	//TestGenAlpha();
	QueryPerformanceCounter(&t3);
	std::cerr << (t3.QuadPart - t2.QuadPart) / tc.QuadPart << "s" << std::endl;
	return 0;
}
