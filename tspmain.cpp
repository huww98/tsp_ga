/*
The GAlib-based genetic algorithm code for the Travelling Salesman Problem (TSP) Berlin52.
作者：wying
单位：华南理工大学软件学院
*/

#include <math.h>
#include <ctime>
#include <string>
#include <iostream>
#include <sstream>
#include "ga/GASStateGA.h"
#include "ga/GASimpleGA.h"
#include "ga/GA1DArrayGenome.h"
#include "ga/garandom.h"

using namespace std;

// Set this up for your favorite TSP.  The sample one is a contrived problem
// with the towns laid out in a grid (so it is easy to figure out what the 
// shortest distance is, and there are many different paths with the same
// shortest path).  File format is that used by the TSPLIB problems.  You can 
// grab more problems from TSPLIB.
// 
constexpr int MAX_TOWNS = 64;
const string TSP_FILE = "berlin52.txt";

int x[MAX_TOWNS], y[MAX_TOWNS];//每个城市的x坐标和y坐标
int DISTANCE[MAX_TOWNS][MAX_TOWNS];//每两个城市之间的旅行成本，是对称的


float TSPObjective(GAGenome&);//计算染色体的旅行总费用的目标函数
void  TSPInitializer(GAGenome&);//TSP问题的染色体初始化算子
int   TSPMutator(GAGenome&, float);//针对TSP问题的染色体变异算子
int   TSPCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);//针对TSP问题的染色体交叉算子

void writeTSPPath(ostream & os, GAGenome& g);//将指定染色体的旅行路线输出到指定文件

GABoolean TSPTerminate(GAGeneticAlgorithm & ga) {
	if (ga.convergence() == 0) return gaFalse;
	return ga.convergence() <= ga.pConvergence() || ga.generation() >= ga.nGenerations() ? gaTrue : gaFalse;
}

int main(int argc, char ** argv) {
	cout << "The GAlib program for the Travelling Salesman Problem (TSP) Berlin52.\n" << endl;

	//从Berlin52.txt文件读出各城市坐标
	double CityID;
	ifstream in(TSP_FILE);
	if (!in) {
		cerr << "could not read data file " << TSP_FILE << "\n";
		exit(1);
	}
	int ntowns = 0;
	do {
		in >> CityID;
		in >> x[ntowns];
		in >> y[ntowns];
		ntowns++;
	} while (!in.eof() && ntowns < MAX_TOWNS);
	in.close();
	if (ntowns >= MAX_TOWNS) {
		cerr << "data file contains more towns than allowed for in the fixed\n";
		cerr << "arrays.  Recompile the program with larger arrays or try a\n";
		cerr << "smaller problem.\n";
		exit(1);
	}

	//计算任意两个城市间的旅行成本
	double dx, dy;
	for (int i = 0; i < ntowns; i++) {
		for (int j = i; j < ntowns; j++) {
			dx = x[i] - x[j]; dy = y[i] - y[j];
			DISTANCE[j][i] = DISTANCE[i][j] = floor(0.5 + sqrt(dx*dx + dy * dy));//注意取四舍五入之后的整数值
		}
	}

	//定义TSP问题的编码方案为一维的整数数组，其固定长度为城市个数
	GA1DArrayGenome<int> genome(ntowns);

	genome.evaluator(::TSPObjective);//为染色体指定计算目标值的函数
	genome.initializer(::TSPInitializer);//为染色体指定自定义的初始化算子
	genome.crossover(::TSPCrossover);//为染色体指定自定义的交叉算子
	genome.mutator(::TSPMutator);//为染色体指定自定义的变异算子

	GASteadyStateGA ga(genome); ga.nReplacement(2); ga.nGenerations(500000); //选用稳态遗传算法进行TSP问题求解，指定其染色体编码方式、每一代要替换的个体数=2、总的运行代数500000，那么搜索的总个体数=2*500000=1000000
	//GASimpleGA ga(genome); ga.elitist(gaTrue); ga.nGenerations(50000);//选用简单遗传算法进行TSP问题求解，采用精英保留策略，指定其染色体编码方式、总的运行代数10000，那么搜索的总个体数=200（种群大小）*5000=1000000
	ga.minimize();//为遗传算法指定优化目的是将目标函数值最小化
	ga.populationSize(2000);//为遗传算法指定种群大小为200
	ga.pMutation(0.02);//为遗传算法指定变异概率
	ga.pCrossover(0.8);//为遗传算法指定交叉概率
	ga.terminator(::TSPTerminate);
	ga.pConvergence(1.0);
	ga.nConvergence(50000);

	unsigned int seed = time(nullptr);
	if (argc > 1)
	{
		stringstream(argv[1]) >> seed;
	}
	cout << "initializing... seed: " << seed << endl;
	ga.initialize(seed);//使用从时钟得到的随机种子初始化遗传算法

	cout << "evolving..." << endl;
	stringstream curveFileName;
	curveFileName << "tspgacurve_" << seed << ".txt";
	std::fstream fgacurve;
	fgacurve.open(curveFileName.str(), std::ios::out);

	//遗传算法开始迭代进化，直到达到指定的代数
	while (!ga.done()) {
		ga.step();//进化一代
		if (ga.generation() % (ga.nGenerations() / 200) == 0)
		{//进化过程中取100个采样点，记录进化过程中的最优目标值收敛信息到文件
			int bestscore = ga.statistics().bestIndividual().score();
			cout << ga.generation() << "\t" << bestscore << endl;
			fgacurve << ga.generation() << "\t" << bestscore << "\n";
		}
	}
	fgacurve.close();

	//遗传算法迭代终止后输出找到的最优旅行路线到文件
	auto & statistics = ga.statistics();
	cout << statistics << endl;
	genome = statistics.bestIndividual();
	//cout << "\n" << "the shortest path found is "  << "\n";
	//writeTSPPath(cout, genome);
	std::fstream fbestpath;
	stringstream pathFileName;
	pathFileName << "tsppath_" << seed << ".txt";
	fbestpath.open(pathFileName.str(), std::ios::out);
	writeTSPPath(fbestpath, genome);
	fbestpath.close();
	cout << "the distance of the shortest path found: " << genome.score() << endl;

	return 0;
}


// Here are the genome operators that we want to use for this problem.
//计算染色体的旅行总费用的目标函数
float TSPObjective(GAGenome& g) {
	GA1DArrayGenome<int> & genome = (GA1DArrayGenome<int> &)g;
	int genomelength = genome.size();//genome.size()获取染色体的长度
	float dist = 0;
	int xx;
	int yy;

	for (int i = 0; i < genomelength; i++) {
		xx = genome.gene(i);
		yy = genome.gene((i + 1) % genomelength);
		dist += DISTANCE[xx - 1][yy - 1];
	}

	return dist;
}

//TSP问题的染色体初始化算子
void TSPInitializer(GAGenome& g) {
	GA1DArrayGenome<int> &genome = dynamic_cast<GA1DArrayGenome<int> &>(g);

	int genomelength = genome.size();
	int i, town;
	bool visit[MAX_TOWNS]{ false };

	for (i = 0; i < genomelength; i++) {
		do {
			town = GARandomInt(1, genomelength);//GARandomInt(1,genomelength)生成1到genomelength之间的一个均匀随机整数
		} while (visit[town - 1]);
		visit[town - 1] = true;
		genome.gene(i, town);//genome.gene(0, town)设置该染色体第0个基因位上的基因值为town
	}
}

//针对TSP问题的染色体变异算子，pmut为变异概率
int TSPMutator(GAGenome& g, float pmut) {
	GA1DArrayGenome<int> &genome = dynamic_cast<GA1DArrayGenome<int> &>(g);

	int genomelength = genome.size();
	float nmutator = 0;//要改变的边的数量
	for (size_t i = 0; i < genomelength; i++)
	{
		if (GARandomFloat() < pmut)
		{
			nmutator++;
		}
	}
	int imutator = 0;
	while (nmutator - imutator >= 2) {
		if (GARandomFloat() < 0.5) {//GARandomFloat()生成0到1之间的一个均匀随机浮点数
		  //以0.5概率使用相互交换变异
			int swapIndex1 = GARandomInt(0, genomelength - 1);
			int swapIndex2 = GARandomInt(0, genomelength - 1);
			int tmp;
			tmp = genome.gene(swapIndex2);
			genome.gene(swapIndex2, genome.gene(swapIndex1));
			genome.gene(swapIndex1, tmp);// swap only one time
			imutator += 4;
		}
		else
		{
			//以0.5概率使用反转变异
			unsigned int inversion_start = GARandomInt(0, genomelength - 1);
			unsigned int inversion_end = GARandomInt(0, genomelength - 1);
			if (inversion_start > inversion_end)
			{
				swap(inversion_start, inversion_end);
			}

			for (; inversion_start < inversion_end; inversion_start++, inversion_end--)
			{
				genome.swap(inversion_start, inversion_end);
			}
			imutator += 2;
		}
	}

	return imutator;
}


void HalfCrossover(const GA1DArrayGenome<int> & parent1,
	const GA1DArrayGenome<int> & parent2,
	GA1DArrayGenome<int> & child,
	unsigned int a, unsigned int b) 
{
	auto l = b - a;
	child.copy(parent2, a, a, l);
	bool visited[MAX_TOWNS]{ false };
	for (size_t i = a; i < b; i++)
	{
		visited[child.gene(i)] = true;
	}
	int childi = 0;
	for (size_t i = 0; i < parent1.size(); i++)
	{
		if (childi == a)
		{
			childi = b;
		}

		auto parentGene = parent1.gene(i);
		if (!visited[parentGene])
		{
			visited[parentGene] = true;
			child.gene(childi, parentGene);
			childi++;
		}
	}
}

//针对TSP问题的染色体交叉算子
int TSPCrossover(const GAGenome& g1, const GAGenome& g2, GAGenome* c1, GAGenome* c2) {
	auto & parent1 = dynamic_cast<const GA1DArrayGenome<int> &>(g1);
	auto & parent2 = (GA1DArrayGenome<int> &)g2;

	int genomelength = parent1.size();

	int nc = 0;

	unsigned int a = GARandomInt(0, parent1.size());
	unsigned int b = GARandomInt(0, parent2.size());
	if (a > b)
	{
		swap(a, b);
	}
	auto l = b - a;

	if (c1) {
		auto & child1 = dynamic_cast<GA1DArrayGenome<int> &>(*c1);
		HalfCrossover(parent1, parent2, child1, a, b);
		nc++;
	}
	if (c2) { 
		auto & child2 = dynamic_cast<GA1DArrayGenome<int> &>(*c2);
		HalfCrossover(parent1, parent2, child2, a, b);
		nc++;
	}

	return nc;
}

//将指定染色体的旅行路线输出到指定文件
void writeTSPPath(ostream & os, GAGenome& g)
{
	GA1DArrayGenome<int> & genome = (GA1DArrayGenome<int> &)g;
	int genomelength = genome.size();
	for (int i = 0; i < genomelength; i++)
	{
		int xx = genome.gene(i);
		os << xx << "    " << x[xx - 1] << "      " << y[xx - 1] << "\n";
	}
}
