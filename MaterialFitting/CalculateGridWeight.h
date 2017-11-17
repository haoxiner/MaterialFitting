#pragma once

#include <iostream>
#include <random>
#include <fstream>
#include <Windows.h>
#include <omp.h>
#include "Mathematics.h"
#undef max

constexpr int WIDTH = 4096;
constexpr int NUM_OF_GRID = WIDTH * WIDTH;
constexpr int SAMPLES_PER_GRID = 1024;
const double SQRT2 = std::sqrt(2.0);
const double SIN_45_DEGREE = SQRT2 / 2.0;

const int NUM_OF_THREAD = 7;

int ClampInsideGrid(int value)
{
	if (value < 0) {
		return 0;
	}
	if (value >= WIDTH) {
		return WIDTH - 1;
	}
	return value;
}

void CalculateGridWeight()
{
	std::ofstream gridSampleFile(R"(F:/gridSample.bin)", std::ios::binary);
	std::ofstream gridVectorFile(R"(F:/gridVector.bin)", std::ios::binary);

	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
	std::cerr << "Running..." << std::endl;

	std::vector<unsigned int> grid(NUM_OF_GRID);
	std::vector<hx::Float3> gridVector(NUM_OF_GRID);

	omp_set_num_threads(NUM_OF_THREAD);

	std::cerr << "NUM OF THREAD: " << omp_get_num_threads() << std::endl;
	std::vector<std::vector<unsigned int>> subGrid(NUM_OF_THREAD);
	for (auto&& v : subGrid) {
		v.resize(NUM_OF_GRID);
	}
	std::vector<std::vector<hx::Float3>> subGridVector(NUM_OF_THREAD);
	for (auto&& v : subGridVector) {
		v.resize(NUM_OF_GRID);
	}
	std::vector<std::default_random_engine> generator(NUM_OF_THREAD);

#pragma omp parallel for
	for (int j = 0; j < SAMPLES_PER_GRID; j++) {
		auto threadID = omp_get_thread_num();
		std::uniform_real_distribution<hx::Float> distributionTheta(0.0, std::nextafter(1.0, std::numeric_limits<hx::Float>::max()));
		std::uniform_real_distribution<hx::Float> distributionPhi(0.0, 1.0);
		for (int i = 0; i < NUM_OF_GRID; i++) {
			auto x = distributionTheta(generator[threadID]);
			auto y = distributionPhi(generator[threadID]);
			auto p3 = hx::UniformSampleSphere(x, y);
			if (p3.z_ < 0.0) {
				std::cerr << "not upper sphere!!" << std::endl;
			}
			auto p2 = hx::MapUpperSphereToPlane(p3);
			hx::Float2 pGrid = { SIN_45_DEGREE * (p2.x_ - p2.y_) + SIN_45_DEGREE, SIN_45_DEGREE * (p2.x_ + p2.y_) + SIN_45_DEGREE };
			int px = ClampInsideGrid(static_cast<int>(pGrid.x_ / SQRT2 * WIDTH));
			int py = ClampInsideGrid(static_cast<int>(pGrid.y_ / SQRT2 * WIDTH));
			subGrid[threadID][px + py * WIDTH] += 1;
			subGridVector[threadID][px + py * WIDTH] += p3;
			if (subGrid[threadID][px + py * WIDTH] == std::numeric_limits<unsigned int>::max()) {
				std::cerr << "Over Flow ERROR !!!!" << std::endl;
			}
		}
	}
	for (int i = 0; i < NUM_OF_THREAD; i++) {
		for (int j = 0; j < NUM_OF_GRID; j++) {
			auto tmp = grid[j] + subGrid[i][j];
			if (tmp < grid[j]) {
				std::cerr << "Over Flow 22222" << std::endl;
			}
			grid[j] = tmp;
			gridVector[j] += subGridVector[i][j];
		}
	}
	for (int i = 0; i < NUM_OF_GRID; i++) {
		if (grid[i] == 0) {
			std::cerr << "Grid has no samples!!" << std::endl;
		}
		gridVector[i] *= (1.0 / grid[i]);
		gridVector[i] = hx::Normalize(gridVector[i]);
	}
	QueryPerformanceCounter(&t2);
	std::cerr << (t2.QuadPart - t1.QuadPart) / tc.QuadPart << "s" << std::endl;

	gridSampleFile.write((const char*)grid.data(), sizeof(grid[0]) * grid.size());
	gridSampleFile.close();
	gridVectorFile.write((const char*)gridVector.data(), sizeof(gridVector[0]) * gridVector.size());
	gridVectorFile.close();
}
