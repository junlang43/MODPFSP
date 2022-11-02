#ifndef __INSGA2_H_
#define __INSGA2_H_

#include "individual.h"


class TINSGA
{
public:
	TINSGA();
	virtual ~TINSGA();

	void init_population();
	void makenewpop();
	TIndividual tournament(TIndividual* pop);
	void fast_nondominated_sort();
	void calu_crowding_distance(int i);
	void crowding_sort(int i);
	void localSearch();

	void run();

	void showPt();
	void showQt();
	void showPA();

	TIndividual* Pt;
	TIndividual* Qt;
	TIndividual** Ft = new TIndividual * [2 * popsize];

	vector<TIndividual>PA;

	int Pnum;
	int Qnum;

	int len[2 * popsize];
	int len_f;
};

TINSGA::TINSGA()
{
	Pnum = 0;
	Qnum = 0;
	Pt = new TIndividual[popsize];
	Qt = new TIndividual[2 * popsize];
	for (int i = 0; i < popsize; i++)
		Ft[i] = new TIndividual[2 * popsize];
	len_f = 0;
	memset(len, 0, sizeof(len));
}

TINSGA::~TINSGA()
{
	delete[] Pt;
	delete[] Qt;
	for (int i = 0; i < popsize; i++)
		delete[] Ft[i];
}

void TINSGA::run()
{
	clock_t startTime, endTime;
	startTime = clock();
	endTime = startTime;
	init_population();
	while ((endTime - startTime) < 0.5 * jobsize * 1000) {
		makenewpop();
		fast_nondominated_sort();
		localSearch();
		fast_nondominated_sort();
		Pnum = 0;
		int maxRank = 0;
		while (Pnum + len[maxRank] <= popsize) {
			for (size_t j = 0; j < len[maxRank]; j++)
				Pt[Pnum++] = Ft[maxRank][j];
			maxRank++;
			if (maxRank >= len_f)break;
		}
		if (maxRank < len_f && Pnum < popsize)
		{
			calu_crowding_distance(maxRank);
			crowding_sort(maxRank);
		}
		int num = Pnum;
		for (size_t j = 0; j < popsize - num; j++)
			Pt[Pnum++] = Ft[maxRank][j];
		for (size_t i = 0; i < popsize; i++)
			Pt[i].updatePA(PA);
		endTime = clock();
		/*showPt();
		cout << endl;
		showPA();
		cout << endl;*/
	}
}

void TINSGA::init_population()
{
	Pt[0].init_MaxSMinTEC();
	Pt[1].init_MaxSMinTFT();
	Pt[2].init_RandSMinTEC();
	Pt[3].init_RandSMinTFT();
	for (size_t i = 4; i < popsize; ++i)
		Pt[i].init_rand();
	Pnum = popsize;
	Qnum = 0;
}

void TINSGA::makenewpop()
{
	TIndividual* X = new TIndividual[popsize];
	for (int i = 0; i < popsize; i++) {
		X[i] = Pt[i];
		int rm = rand_int(0, 2);
		switch (rm) {
		case 0:
			X[i].insertNew(PA); break;
		case 1:
			X[i].swapNew(PA); break;
		case 2:
			X[i].hybridNew(PA); break;
		}
	}
	for (int i = 0; i < popsize; i++) {
		if (rand_real(0, 1) < 0.5)
			Qt[i] = tournament(Pt);
		else
			Qt[i] = tournament(X);
		int rm = rand_int(0, 2);
		switch (rm) {
		case 0:
			Qt[i].insertNew(PA); break;
		case 1:
			Qt[i].swapNew(PA); break;
		case 2:
			Qt[i].hybridNew(PA); break;
		}
	}
	for (size_t i = 0; i < popsize; i++)
		Qt[i + popsize] = X[i];
}

void TINSGA::fast_nondominated_sort()
{
	int i;
	TIndividual* H = new TIndividual[2 * popsize];
	int h_len = 0;
	len_f = 0;
	for (i = 0; i < 2 * popsize; i++)
	{
		Qt[i].np = 0;
		Qt[i].dominated_num = 0;
		len[i] = 0;
	}
	for (i = 0; i < 2 * popsize; i++)
	{
		for (int j = 0; j < 2 * popsize; j++)
		{
			if (i != j)
			{
				if (is_dominated(Qt[i].y_obj, Qt[j].y_obj))
					Qt[i].sp[Qt[i].dominated_num++] = j;
				else if (is_dominated(Qt[j].y_obj, Qt[i].y_obj))
					Qt[i].np += 1;
			}
		}
		if (Qt[i].np == 0)
		{
			len_f = 1;
			Ft[0][len[0]++] = Qt[i];
		}

	}
	i = 0;
	while (len[i] != 0)
	{
		h_len = 0;
		for (int j = 0; j < len[i]; j++)
		{
			for (int k = 0; k < Ft[i][j].dominated_num; k++)
			{
				Qt[Ft[i][j].sp[k]].np--;
				if (Qt[Ft[i][j].sp[k]].np == 0)
				{
					H[h_len++] = Qt[Ft[i][j].sp[k]];
					Qt[Ft[i][j].sp[k]].rank = i + 2;
				}
			}
		}
		i++;
		len[i] = h_len;
		if (h_len != 0)
		{
			len_f++;
			for (int j = 0; j < len[i]; j++) {
				Ft[i][j] = H[j];
			}
		}
	}
}

void TINSGA::calu_crowding_distance(int i)
{
	int n = len[i];
	double m_max, m_min;
	int j;
	for (j = 0; j < n; j++)
		Ft[i][j].crowding_distance = 0;
	Ft[i][0].crowding_distance = Ft[i][n - 1].crowding_distance = 0xffffff;
	qsort(Ft[i], n, sizeof(TIndividual), cmp1);
	m_max = -0xfffff;
	m_min = 0xfffff;
	if (m_max < Ft[i][n - 1].y_obj[0])
		m_max = Ft[i][n - 1].y_obj[0];
	if (m_min > Ft[i][0].y_obj[0])
		m_min = Ft[i][0].y_obj[0];
	for (j = 1; j < n - 1; j++) {
		Ft[i][j].crowding_distance += (Ft[i][j + 1].y_obj[0] - Ft[i][j - 1].y_obj[0]) / (m_max - m_min);
	}

	qsort(Ft[i], n, sizeof(TIndividual), cmp2);
	m_max = -0xfffff;
	m_min = 0xfffff;
	if (m_max < Ft[i][n - 1].y_obj[1])
		m_max = Ft[i][n - 1].y_obj[1];
	if (m_min > Ft[i][0].y_obj[1])
		m_min = Ft[i][0].y_obj[1];
	for (j = 1; j < n - 1; j++)
		Ft[i][j].crowding_distance += (Ft[i][j + 1].y_obj[1] - Ft[i][j - 1].y_obj[1]) / (m_max - m_min);
}
void TINSGA::crowding_sort(int i)
{
	int n;
	n = len[i];
	qsort(Ft[i], n, sizeof(TIndividual), cmp_d);
}

void TINSGA::localSearch()
{
	Qnum = 0;
	for (size_t i = 0; i < len_f; i++) {
		for (size_t j = 0; j < len[i]; j++)
			Qt[Qnum++] = Ft[i][j];
	}
	for (size_t i = 0; i < len[0]; i++) {
 		Qt[i].localIntensifization(PA);
	}
}

TIndividual TINSGA::tournament(TIndividual* pop)
{
	int p1 = rand_int(0, popsize - 1);
	int p2 = p1;
	while (p2 == p1)
		p2 = rand_int(0, popsize - 1);
	if (is_dominated(pop[p1].y_obj, pop[p2].y_obj))
		return pop[p1];
	else if (is_dominated(pop[p2].y_obj, pop[p1].y_obj))
		return pop[p2];
	else
		return rand_real(0, 1) < 0.5 ? pop[p1] : pop[p2];
}

void TINSGA::showPt()
{
	for (int i = 0; i < popsize; i++) {
		Pt[i].show_objective();
		Pt[i].show_variable();
	}
}

void TINSGA::showQt()
{
	for (int i = 0; i < popsize; i++) {
		Qt[i].show_objective();
		Qt[i].show_variable();
	}
}

void TINSGA::showPA()
{
	for (auto it = PA.begin(); it != PA.end();it++) {
		it->show_objective();
		//it->show_variable();
	}
}

#endif

