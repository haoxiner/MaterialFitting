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
	// grid weight
	std::ifstream gridWeightFile(R"(F:/gridWeight.bin)", std::ios::binary);
	std::vector<hx::Float> gridWeight(NUM_OF_GRID);
	gridWeightFile.read((char*)gridWeight.data(), gridWeight.size() * sizeof(gridWeight[0]));
	gridWeightFile.close();

	// grid center
	//std::vector<hx::Float3> gridCenter(NUM_OF_GRID);
	//std::ifstream gridCenterFile(R"(F:/gridCenter.bin)", std::ios::binary);
	//gridCenterFile.read((char*)gridCenter.data(), gridCenter.size() * sizeof(gridCenter[0]));
	//gridCenterFile.close();

	// grid ward pdf
	std::vector<hx::Float> gridWardPdf(NUM_OF_GRID);
	
	std::vector<hx::Float> alphas = { 0.5 };
	for (auto alpha : alphas) {
		hx::MapPlaneToUpperSphere({});
		for (int degree = 0; degree < 90; degree++) {
			hx::Float cosTheta = std::cos(hx::DegreeToRadian(static_cast<hx::Float>(degree)));
			hx::Float3 wo = { std::sqrt(1.0 - cosTheta * cosTheta), 0.0, cosTheta };

			std::vector<hx::Float> pdfX(WIDTH);
			std::vector<hx::Float> pdfYAfterX(WIDTH);
			for (int y = 0; y < WIDTH; y++) {
				for (int x = 0; x < WIDTH; x++) {
					int idx = x + y * WIDTH;
					hx::Float px = (x + 0.5) / WIDTH * SQRT2;
					hx::Float py = (y + 0.5) / WIDTH * SQRT2;
					auto gridCenter = hx::MapPlaneToUpperSphere({ SIN_45_DEGREE * (px + py) - 1.0, SIN_45_DEGREE * (py - px) });
					gridWardPdf[idx] = hx::Ward(alpha, gridCenter, wo) * gridWeight[idx];

					pdfX[x] += gridWardPdf[idx];
					pdfYAfterX[y] += gridWardPdf[idx];
				}
			}
			std::vector<hx::Float> cdfX(WIDTH);
			std::vector<hx::Float> cdfYAfterX(WIDTH);
			cdfX[0] = pdfX[0];
			cdfYAfterX[0] = pdfYAfterX[0];
			for (int i = 1; i < WIDTH; i++) {
				cdfX[i] = cdfX[i - 1] + pdfX[i];
				cdfYAfterX[i] = cdfYAfterX[i - 1] + pdfYAfterX[i];
			}
			std::cerr << std::accumulate(gridWeight.begin(), gridWeight.end(), 0.0) << std::endl;
		}
	}
}

int main()
{
	/*CalculateGridWeight();
	GenGridPdf();*/
	//GenGridCenter();
	GenGridWardPdf();
	return 0;
}
