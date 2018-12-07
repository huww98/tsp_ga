/*
The GAlib-based genetic algorithm code for the Travelling Salesman Problem (TSP) Berlin52.
���ߣ�wying
��λ����������ѧ���ѧԺ
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

int x[MAX_TOWNS], y[MAX_TOWNS];//ÿ�����е�x�����y����
int DISTANCE[MAX_TOWNS][MAX_TOWNS];//ÿ��������֮������гɱ����ǶԳƵ�


float TSPObjective(GAGenome&);//����Ⱦɫ��������ܷ��õ�Ŀ�꺯��
void  TSPInitializer(GAGenome&);//TSP�����Ⱦɫ���ʼ������
int   TSPMutator(GAGenome&, float);//���TSP�����Ⱦɫ���������
int   TSPCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);//���TSP�����Ⱦɫ�彻������

void writeTSPPath(ostream & os, GAGenome& g);//��ָ��Ⱦɫ�������·�������ָ���ļ�

GABoolean TSPTerminate(GAGeneticAlgorithm & ga) {
	if (ga.convergence() == 0) return gaFalse;
	return ga.convergence() <= ga.pConvergence() || ga.generation() >= ga.nGenerations() ? gaTrue : gaFalse;
}

int main(int argc, char ** argv) {
	cout << "The GAlib program for the Travelling Salesman Problem (TSP) Berlin52.\n" << endl;

	//��Berlin52.txt�ļ���������������
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

	//���������������м�����гɱ�
	double dx, dy;
	for (int i = 0; i < ntowns; i++) {
		for (int j = i; j < ntowns; j++) {
			dx = x[i] - x[j]; dy = y[i] - y[j];
			DISTANCE[j][i] = DISTANCE[i][j] = floor(0.5 + sqrt(dx*dx + dy * dy));//ע��ȡ��������֮�������ֵ
		}
	}

	//����TSP����ı��뷽��Ϊһά���������飬��̶�����Ϊ���и���
	GA1DArrayGenome<int> genome(ntowns);

	genome.evaluator(::TSPObjective);//ΪȾɫ��ָ������Ŀ��ֵ�ĺ���
	genome.initializer(::TSPInitializer);//ΪȾɫ��ָ���Զ���ĳ�ʼ������
	genome.crossover(::TSPCrossover);//ΪȾɫ��ָ���Զ���Ľ�������
	genome.mutator(::TSPMutator);//ΪȾɫ��ָ���Զ���ı�������

	GASteadyStateGA ga(genome); ga.nReplacement(2); ga.nGenerations(500000); //ѡ����̬�Ŵ��㷨����TSP������⣬ָ����Ⱦɫ����뷽ʽ��ÿһ��Ҫ�滻�ĸ�����=2���ܵ����д���500000����ô�������ܸ�����=2*500000=1000000
	//GASimpleGA ga(genome); ga.elitist(gaTrue); ga.nGenerations(50000);//ѡ�ü��Ŵ��㷨����TSP������⣬���þ�Ӣ�������ԣ�ָ����Ⱦɫ����뷽ʽ���ܵ����д���10000����ô�������ܸ�����=200����Ⱥ��С��*5000=1000000
	ga.minimize();//Ϊ�Ŵ��㷨ָ���Ż�Ŀ���ǽ�Ŀ�꺯��ֵ��С��
	ga.populationSize(2000);//Ϊ�Ŵ��㷨ָ����Ⱥ��СΪ200
	ga.pMutation(0.02);//Ϊ�Ŵ��㷨ָ���������
	ga.pCrossover(0.8);//Ϊ�Ŵ��㷨ָ���������
	ga.terminator(::TSPTerminate);
	ga.pConvergence(1.0);
	ga.nConvergence(50000);

	unsigned int seed = time(nullptr);
	if (argc > 1)
	{
		stringstream(argv[1]) >> seed;
	}
	cout << "initializing... seed: " << seed << endl;
	ga.initialize(seed);//ʹ�ô�ʱ�ӵõ���������ӳ�ʼ���Ŵ��㷨

	cout << "evolving..." << endl;
	stringstream curveFileName;
	curveFileName << "tspgacurve_" << seed << ".txt";
	std::fstream fgacurve;
	fgacurve.open(curveFileName.str(), std::ios::out);

	//�Ŵ��㷨��ʼ����������ֱ���ﵽָ���Ĵ���
	while (!ga.done()) {
		ga.step();//����һ��
		if (ga.generation() % (ga.nGenerations() / 200) == 0)
		{//����������ȡ100�������㣬��¼���������е�����Ŀ��ֵ������Ϣ���ļ�
			int bestscore = ga.statistics().bestIndividual().score();
			cout << ga.generation() << "\t" << bestscore << endl;
			fgacurve << ga.generation() << "\t" << bestscore << "\n";
		}
	}
	fgacurve.close();

	//�Ŵ��㷨������ֹ������ҵ�����������·�ߵ��ļ�
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
//����Ⱦɫ��������ܷ��õ�Ŀ�꺯��
float TSPObjective(GAGenome& g) {
	GA1DArrayGenome<int> & genome = (GA1DArrayGenome<int> &)g;
	int genomelength = genome.size();//genome.size()��ȡȾɫ��ĳ���
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

//TSP�����Ⱦɫ���ʼ������
void TSPInitializer(GAGenome& g) {
	GA1DArrayGenome<int> &genome = dynamic_cast<GA1DArrayGenome<int> &>(g);

	int genomelength = genome.size();
	int i, town;
	bool visit[MAX_TOWNS]{ false };

	for (i = 0; i < genomelength; i++) {
		do {
			town = GARandomInt(1, genomelength);//GARandomInt(1,genomelength)����1��genomelength֮���һ�������������
		} while (visit[town - 1]);
		visit[town - 1] = true;
		genome.gene(i, town);//genome.gene(0, town)���ø�Ⱦɫ���0������λ�ϵĻ���ֵΪtown
	}
}

//���TSP�����Ⱦɫ��������ӣ�pmutΪ�������
int TSPMutator(GAGenome& g, float pmut) {
	GA1DArrayGenome<int> &genome = dynamic_cast<GA1DArrayGenome<int> &>(g);

	int genomelength = genome.size();
	float nmutator = 0;//Ҫ�ı�ıߵ�����
	for (size_t i = 0; i < genomelength; i++)
	{
		if (GARandomFloat() < pmut)
		{
			nmutator++;
		}
	}
	int imutator = 0;
	while (nmutator - imutator >= 2) {
		if (GARandomFloat() < 0.5) {//GARandomFloat()����0��1֮���һ���������������
		  //��0.5����ʹ���໥��������
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
			//��0.5����ʹ�÷�ת����
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

//���TSP�����Ⱦɫ�彻������
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

//��ָ��Ⱦɫ�������·�������ָ���ļ�
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
