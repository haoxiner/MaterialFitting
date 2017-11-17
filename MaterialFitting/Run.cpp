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

void GenGridWardPdf()
{
	std::ifstream file(R"(F:/weight.bin)", std::ios::binary);
	std::ofstream out(R"(F:/wardpdf.bin)", std::ios::binary);
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

int main()
{
	/*CalculateGridWeight();
	GenGridPdf();*/
	std::ifstream file(R"(F:/gridVector.bin)", std::ios::binary);
	
	std::vector<hx::Float3> src(NUM_OF_GRID);
	file.read((char*)src.data(), src.size() * sizeof(hx::Float3));
	file.close();

	int cnt = 0;
	for (auto v : src) {
		if (v.x_ < 0.0 || v.y_ < 0.0) {
			cnt++;
			std::cerr << v.x_ << ", " << v.y_ << ", " << v.z_ << std::endl;
		}
		if (cnt > 30) {
			break;
		}
	}
	return 0;
}
