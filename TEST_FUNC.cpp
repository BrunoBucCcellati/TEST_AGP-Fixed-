#include "pch.h"
#include <random>
#define M_PI 3.14159265358979323846
extern "C" __declspec(dllexport)double HillFunc(double x, time_t now)
{
	const unsigned short HH_COUNT = 15;
	double a = 0, b = 1;
	srand((unsigned short)now);
	std::default_random_engine generator(rand());
	std::uniform_real_distribution<> distribution(a, b);
	double res = -1.1 + distribution(generator) * 2;
	unsigned short i = 1;
	while (i < HH_COUNT)
	{
		res += (-1.1 + distribution(generator) * 2) * sin(M_PI * i * x * 2) + (-1.1 + distribution(generator) * 2) * cos(M_PI * i * x * 2);
		i++;
	}
	return res;
}
extern "C" __declspec(dllexport)double ShekelFunc(double x, time_t now)
{
	const unsigned short SH_COUNT = 11;
	double a = 0, b = 1, res = 0;
	srand((unsigned short)now);
	std::default_random_engine generator(rand());
	std::uniform_real_distribution<> distribution(a, b);
	unsigned short i = 1;
	while (i < SH_COUNT)
	{
		res -= 1 / ((5 + 20 * distribution(generator)) * pow(x - 10 * distribution(generator), 2) + 1 + 0.4 * distribution(generator));
		i++;
	}
	return res;
}
extern "C" __declspec(dllexport) double GrishaginFunc(double x1, double x2, time_t now)
{
	const unsigned short GR_COUNT = 8;
	double a = -1, b = 1, part1 = 0, part2 = 0;
	srand((unsigned short)now);
	std::default_random_engine generator(rand());
	std::uniform_real_distribution<> distribution(a, b);
	unsigned short i = 1, j;
	while (i < GR_COUNT)
	{
		j = 1;
		while (j < GR_COUNT)
		{
			part1 += distribution(generator) * sin(M_PI * i * x1) * sin(M_PI * j * x2) + distribution(generator) * cos(M_PI * i * x1) * cos(M_PI * j * x2);
			part2 += distribution(generator) * sin(M_PI * i * x1) * sin(M_PI * j * x2) - distribution(generator) * cos(M_PI * i * x1) * cos(M_PI * j * x2);
			j++;
		}
		i++;
	}
	return -sqrt(part1 * part1 + part2 * part2);
}