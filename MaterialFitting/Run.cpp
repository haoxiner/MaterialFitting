#include "CalculateGridWeight.h"

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

void GenGridWardPdf()
{
	std::ifstream file(R"(F:/gridWeight.bin)", std::ios::binary);
	std::vector<unsigned int> gridWeight(NUM_OF_GRID);
	file.read((char*)gridWeight.data(), gridWeight.size() * sizeof(unsigned int));
	file.close();

	std::vector<hx::Float3> gridCenter(NUM_OF_GRID);
	for (int y = 0; y < WIDTH; y++) {
		for (int x = 0; x < WIDTH; x++) {
			gridCenter[x + y * WIDTH] = hx::MapPlaneToUpperSphere({ SIN_45_DEGREE * (x + y + 1.0) - 1.0, SIN_45_DEGREE * (y - x) });
		}
	}

	std::vector<hx::Float> alphas = { 0.007 };
	for (auto alpha : alphas) {
		hx::MapPlaneToUpperSphere({});
		for (int degree = 0; degree < 90; degree++) {

		}
	}
}

int main()
{
	/*CalculateGridWeight();
	GenGridPdf();*/
	GenGridCenter();
	return 0;
}
