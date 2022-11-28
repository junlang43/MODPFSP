#ifndef __INDIVIDUAL_H_
#define __INDIVIDUAL_H_

#include "global.h"


class TIndividual {
public:
	TIndividual();
	virtual ~TIndividual();

	vector <vector<int>> x_var;//solution
	vector<vector<int>>x_speed;//speed n*m
	vector <double> y_obj;//object
	vector<vector<pair<int, int>>>critial_path;

	double PEC[factorysize];
	double SEC[factorysize];
	double TEC[factorysize];
	double Compelet[factorysize];
	double idleTime[jobsize][machinesize];//n*m for x_speed

	int fmax;
	int ftec;

	int sp[2 * popsize];//被支配个体集合SP。该量是可行解空间中所有被个体p支配的个体组成的集合。
	int dominated_num;//集合sp的个数
	int np;//支配个数np。该量是在可行解空间中可以支配个体p的所有个体的数量。		
	int rank;//优先级，Pareto级别为当前最高级
	double crowding_distance;//拥挤距离

	void factory_obj(int f);
	void obj_eval();
	void speedup(int f);
	void randspeedup(int f);
	void speeddown(int f);
	void randspeeddown(int f);

	void init_rand();
	void init_MaxSMinTFT();
	void init_MaxSMinTEC();
	void init_RandSMinTFT();
	void init_RandSMinTEC();


	void insertNew(vector<vector<double>> &archive);
	void swapNew(vector<vector<double>> &archive);
	void hybridNew(vector<vector<double>> &archive);

	void updatePA(vector<vector<double>> &archive);

	void localIntensifization(vector<vector<double>> &archive);

	void show_objective();
	void show_variable();
	void show_speed();
};

TIndividual::TIndividual()
{
	x_var.resize(factorysize);
	x_speed.resize(jobsize);
	critial_path.resize(factorysize);
	memset(idleTime, 0, sizeof(idleTime));
	for (size_t n = 0; n < numObjectives; ++n)
		y_obj.push_back(0.0);
	for (size_t n = 0; n < factorysize; ++n) {
		PEC[n] = 0.0;
		SEC[n] = 0.0;
		Compelet[n] = 0.0;
		fmax = 0;
		ftec = 0;
	}
}

TIndividual::~TIndividual()
{

}

void TIndividual::obj_eval()
{
	for (int i = 0; i < factorysize; i++) {
		if (x_var[i].size() != 0) {
			factory_obj(i);
		}		
	}
	vector<double>tec;
	for (int i = 0; i < factorysize; i++) {
		tec.push_back(PEC[i] + SEC[i]);
	}
	fmax = max_element(Compelet, Compelet + factorysize) - Compelet;
	ftec = max_element(tec.begin(), tec.end()) - tec.begin();
	y_obj[0] = Compelet[fmax];
	y_obj[1] = accumulate(PEC, PEC + factorysize, 0.0) + accumulate(SEC, SEC + factorysize, 0.0);
	
}
void TIndividual::factory_obj(int f)
{
	size_t subjSize = x_var[f].size();
	vector<vector<double>>C(subjSize, vector<double>(machinesize, 0.0));
	int ii = 0;
	double NewT[jobsize][machinesize] = { 0 };
	double sumP[machinesize] = { 0 };
	if (subjSize < 1) {
		cout << "wrong" << endl;
	}
	PEC[f] = 0;
	SEC[f] = 0;
	for (vector<int>::iterator it = x_var[f].begin(); it != x_var[f].end(); it++) {
		size_t n = *it;
		if (n <= 0) {
			cout << "job index wrong" << endl;
			for (auto t = x_var[f].begin(); t != x_var[f].end(); t++)
				cout << *t << endl;
			cout << endl;
		}
		for (size_t j = 0; j < machinesize; j++) {
			NewT[ii][j] = double(Time[n - 1][j]) / Speed[x_speed[n - 1][j]];
			PEC[f] += (NewT[ii][j] * 4 * Speed[x_speed[n - 1][j]] * Speed[x_speed[n - 1][j]]);
			sumP[j] += NewT[ii][j];
		}
		ii++;
	}
	C[0][0] = NewT[0][0];
	for (size_t i = 1; i < subjSize; i++) {
		C[i][0] = C[i - 1][0] + NewT[i][0];
	}
	for (size_t i = 1; i < machinesize; i++) {
		C[0][i] = C[0][i - 1] + NewT[0][i];
	}
	for (size_t i = 1; i < subjSize; i++) {
		for (size_t j = 1; j < machinesize; j++) {
			C[i][j] = max(C[i][j - 1], C[i - 1][j]) + NewT[i][j];
			idleTime[x_var[f][i - 1] - 1][j] = max(C[i][j - 1] - C[i - 1][j], 0);
		}
	}
	for (size_t i = 0; i < machinesize; i++) {
		SEC[f] += (C.back()[i] -C.front()[i]- sumP[i]) * SP;
	}
	Compelet[f] = C.back().back();
	TEC[f] = PEC[f] + SEC[f];
	size_t i = subjSize - 1;
	size_t j = machinesize - 1;
	critial_path[f].clear();
	critial_path[f].push_back(pair<int, int>{x_var[f][i], machinesize - 1});
	while (i >= 0) {
		if (i > 0 && j > 0) {
			if (C[i][j - 1] > C[i - 1][j])
				--j;
			else
				--i;
			critial_path[f].insert(critial_path[f].begin(), pair<int, int>{x_var[f][i], j});
		}
		else if (i == 0 && j > 0) {
			--j;
			critial_path[f].insert(critial_path[f].begin(), pair<int, int>{x_var[f][i], j});
		}
		else if (j == 0 && i > 0) {
			--i;
			critial_path[f].insert(critial_path[f].begin(), pair<int, int>{x_var[f][i], j});
		}
		else {
			break;
		}
	}
}

void TIndividual::speedup(int f)
{
	double cost = TEC[f];
	for (auto it = critial_path[f].begin(); it != critial_path[f].end(); it++) {
		if (x_speed[it->first - 1][it->second] < 4) {
			x_speed[it->first - 1][it->second]++;
		}
	}
}

void TIndividual::randspeedup(int f) {
	for (auto it = critial_path[f].begin(); it != critial_path[f].end(); it++) {
		x_speed[it->first - 1][it->second] += rand_int(0, 4 - x_speed[it->first - 1][it->second]);
	}
}

void TIndividual::speeddown(int f)
{
	double cost = y_obj[0];
	for (auto it = x_var[f].begin(); it != x_var[f].end()-1; it++) {
		size_t n = *it;
		if (n <= 0) {
			cout << "job index wrong" << endl;
			for (auto t = x_var[f].begin(); t != x_var[f].end(); t++)
				cout << *t << endl;
			cout << endl;
		}
		for (int i = 1; i < machinesize; ++i) {
			if (find(critial_path[f].begin(), critial_path[f].end(), pair<int, int>{n, i}) != critial_path[f].end())
				continue;
			while (x_speed[n - 1][i] > 0 ) {
				--x_speed[n - 1][i];
				factory_obj(f);
				y_obj[0] = *max_element(Compelet, Compelet + factorysize);
				if (y_obj[0] > cost) {
					++x_speed[n - 1][i];
					break;
				}
			}
		}
	}
}

void TIndividual::randspeeddown(int f)
{
	for (auto it = x_var[f].begin(); it != x_var[f].end(); it++) {
		size_t n = *it;
		for (int i = 0; i < machinesize; ++i) {
			if (find(critial_path[f].begin(), critial_path[f].end(), pair<int, int>{n, i}) != critial_path[f].end())
				continue;
			x_speed[n - 1][i] -= rand_int(0, x_speed[n - 1][i] - 0);
			if (x_speed[n - 1][i] < 0)	cout << "x_speed wrong" << endl;
		}
	}
}

void  TIndividual::init_rand()
{
	int c = 0;
	vector<size_t>pi(jobsize,0);
	x_speed.resize(jobsize);
	for (size_t i = 0; i < jobsize; ++i) {
		x_speed[i].resize(machinesize);
		pi[i] = i + 1;
		for (size_t j = 0; j < machinesize; ++j) {
			x_speed[i][j] = rand_int(0, 4);
		}
	}
	random_shuffle(pi.begin(), pi.end());
	for (auto it = pi.begin(); it != pi.end(); it++) {
		x_var[c % factorysize].push_back(*it);
		c++;
	}
	obj_eval();
}

void TIndividual::init_MaxSMinTFT()
{
	vector<pair<int, double>> sumProcess;
	int c = 0;
	for (auto it = Time.begin(); it != Time.end(); it++) {
		x_speed[c].resize(machinesize, 4);
		double sum = (double)accumulate((*it).begin(), (*it).end(), 0);
		sum /= Speed[4];
		sumProcess.push_back({ ++c,sum});
	}
	sort(sumProcess.begin(), sumProcess.end(), cmp);
	double cost = 0;
	for (auto it = sumProcess.begin(); it != sumProcess.end(); it++) {
		double complete = DBL_MAX;
		double tempcomplete;
		int factory = 0;
		for (int k = 0; k < factorysize; k++) {
			x_var[k].push_back(it->first);
			factory_obj(k);
			tempcomplete = max(cost,Compelet[k]);
			if (tempcomplete < complete) {
				complete = tempcomplete;
				factory = k;
			}
			x_var[k].pop_back();
		}
		x_var[factory].push_back(it->first);
		cost = complete;
	}
	obj_eval();
	for (int i = 0; i < factorysize; i++) {
		speeddown(i);
	}
	obj_eval();
}

void TIndividual::init_MaxSMinTEC()
{
	vector<pair<int, double>> sumProcess;
	int c = 0;
	for (auto it = Time.begin(); it != Time.end(); it++) {
		x_speed[c].resize(machinesize, 0);
		double sum = (double)accumulate((*it).begin(), (*it).end(), 0);
		sum /= Speed[0];
		sumProcess.push_back({ ++c,sum });
	}
	sort(sumProcess.begin(), sumProcess.end(), cmp);
	for (auto it = sumProcess.begin(); it != sumProcess.end(); it++) {
		double complete = DBL_MAX;
		double deltatec;
		int factory = 0;
		for (int k = 0; k < factorysize; k++) {
			double tec = PEC[k] + SEC[k];
			x_var[k].push_back(it->first);
			factory_obj(k);
			deltatec = PEC[k] + SEC[k] - tec;
			if (deltatec < complete) {
				complete = deltatec;
				factory = k;
			}
			x_var[k].pop_back();
		}
		x_var[factory].push_back(it->first);
		obj_eval();
	}
	for (int i = 0; i < factorysize; i++) {
		speedup(i);
	}
	obj_eval();
}

void TIndividual::init_RandSMinTFT()
{
	vector<pair<int, double>> sumProcess;
	int c = 0;
	for (auto it = Time.begin(); it != Time.end(); it++) {
		x_speed[c].resize(machinesize);
		double sum = 0.0;
		for (size_t i = 0; i < machinesize; ++i) {
			x_speed[c][i] = rand_int(0, 4);
			sum += (*it)[i] / Speed[x_speed[c][i]];
		}
		sumProcess.push_back({ ++c,sum });
	}
	sort(sumProcess.begin(), sumProcess.end(), cmp);
	double cost = 0;
	for (auto it = sumProcess.begin(); it != sumProcess.end(); it++) {
		double complete = DBL_MAX;
		double tempcomplete;
		int factory = 0;
		for (int k = 0; k < factorysize; k++) {
			x_var[k].push_back(it->first);
			factory_obj(k);
			tempcomplete = max(cost, Compelet[k]);
			if (tempcomplete < complete) {
				complete = tempcomplete;
				factory = k;
			}
			x_var[k].pop_back();
		}
		x_var[factory].push_back(it->first);
		cost = complete;
	}
	obj_eval();
	for (int i = 0; i < factorysize; i++) {
		speeddown(i);
	}
	obj_eval();
}

void TIndividual::init_RandSMinTEC()
{
	vector<pair<int, double>> sumProcess;
	int c = 0;
	for (auto it = Time.begin(); it != Time.end(); it++) {
		double sum = 0.0;
		x_speed[c].resize(machinesize);
		for (size_t i = 0; i < machinesize; i++) {
			x_speed[c][i] = rand_int(0, 4);
			sum += (*it)[i] / Speed[x_speed[c][i]];
		}
		sumProcess.push_back({ ++c,sum });
	}
	sort(sumProcess.begin(), sumProcess.end(), cmp);
	for (auto it = sumProcess.begin(); it != sumProcess.end(); it++) {
		double complete = DBL_MAX;
		double deltatec;
		int factory = 0;
		for (int k = 0; k < factorysize; k++) {
			double tec = PEC[k] + SEC[k];
			x_var[k].push_back(it->first);
			factory_obj(k);
			deltatec = PEC[k] + SEC[k] - tec;
			if (deltatec < complete) {
				complete = deltatec;
				factory = k;
			}
			x_var[k].pop_back();
		}
		x_var[factory].push_back(it->first);
		obj_eval();
	}
	for (int i = 0; i < factorysize; i++) {
		speeddown(i);
	}
	obj_eval();
}

void TIndividual::updatePA(vector<vector<double>> &archive) 
{
	bool flag = true;
	for (vector<vector<double>>::iterator it = archive.begin(); it != archive.end();) {
		if (is_dominated(y_obj, *it))
			it = archive.erase(it);
		else {
			if ((y_obj[0] == (*it)[0] && y_obj[1] == (*it)[1]) || is_dominated(*it, y_obj)) {
				flag = false;
				break;
			}
			it++;
		}
	}
	if (flag) {
		archive.push_back(y_obj);
	}
}

void TIndividual::insertNew(vector<vector<double>> &archive)
{
	int objfactory = 0;
	bool flag = 0;
	int job = 0;
	vector<int>jobtried;
	vector <vector<int>> temp_var(x_var);
	vector<vector<int>>temp_speed(x_speed);
	vector<double>temp_obj(y_obj);
	if (rand_real(0, 1) < 0.5)
		flag = 1;
	if (flag == 0)
		objfactory = fmax;
	else
		objfactory = ftec;
	if (x_var[objfactory].size() < 1)
		cout << "factory size wrong" << endl;
	while (jobtried.size() < min(2,ceil(x_var[objfactory].size()/2.0))) {
		x_var.assign(temp_var.begin(), temp_var.end());
		x_speed.assign(temp_speed.begin(), temp_speed.end());
		y_obj.assign(temp_obj.begin(), temp_obj.end());
		int pos = rand_int(0, x_var[objfactory].size() - 1);
		while(find(jobtried.begin(),jobtried.end(), x_var[objfactory][pos])!= jobtried.end())
			pos = rand_int(0, x_var[objfactory].size() - 1);
		job = x_var[objfactory][pos];
		x_var[objfactory].erase(x_var[objfactory].begin() + pos);
		jobtried.push_back(job);
		//factory_obj(objfactory);
		if (flag == 0) {
			randspeedup(objfactory);
			speedup(objfactory);
		}
		else {
			randspeeddown(objfactory);
			speeddown(objfactory);
		}
		for (int i = 0; i < factorysize; i++) {
			if (i != objfactory) {
				for (int j = 0; j < x_var[i].size(); j++) {
					x_speed.assign(temp_speed.begin(), temp_speed.end());
					x_var[i].insert(x_var[i].begin() + j, job);
					//factory_obj(i);
					if (flag == 0) {
						randspeedup(i);
						speedup(i);
					}
					else {
						randspeeddown(i);
						speeddown(i);
					}
					obj_eval();
					if (is_dominated(y_obj, temp_obj))
						return;
					else if (!is_dominated(temp_obj, y_obj))
						updatePA(archive);
					x_var[i].erase(x_var[i].begin() + j);
				}
			}
		}
	}
	x_var[factorysize - 1].push_back(job);
}

void TIndividual::swapNew(vector<vector<double>> &archive)
{
	int objfactory = 0;
	bool flag = 0;
	int job = 0;
	int pos = 0;
	vector<int>jobtried;
	vector <vector<int>> temp_var(x_var);
	vector<vector<int>>temp_speed(x_speed);
	vector<double>temp_obj(y_obj);
	if (rand_real(0, 1) < 0.5)
		flag = 1;
	if (flag == 0)
		objfactory = fmax;
	else
		objfactory = ftec;
	if (x_var[objfactory].size() < 1)
		cout << "factory size wrong" << endl;
	while (jobtried.size() < min(2, ceil(x_var[objfactory].size() / 2.0))) {
		x_var.assign(temp_var.begin(), temp_var.end());
		y_obj.assign(temp_obj.begin(), temp_obj.end());
		pos = rand_int(0, x_var[objfactory].size() - 1);
		while (find(jobtried.begin(), jobtried.end(), x_var[objfactory][pos]) != jobtried.end())
			pos = rand_int(0, x_var[objfactory].size() - 1);
		job = x_var[objfactory][pos];
		x_var[objfactory].erase(x_var[objfactory].begin() + pos);
		jobtried.push_back(job);
		for (int i = 0; i < factorysize; i++) {
			if (i != objfactory) {
				for (int j = 0; j < x_var[i].size(); j++) {
					x_var.assign(temp_var.begin(), temp_var.end());
					x_speed.assign(temp_speed.begin(), temp_speed.end());
					x_var[objfactory][pos] = x_var[i][j];
					x_var[i][j] = job;
					//factory_obj(objfactory);
					//factory_obj(i);
					if (flag == 0) {
						randspeedup(objfactory);
						speedup(objfactory);
						randspeedup(i);
						speedup(i);
					}
					else {
						randspeeddown(objfactory);
						speeddown(objfactory);
						randspeeddown(i);
						speeddown(i);
					}
					obj_eval();
					if (is_dominated(y_obj, temp_obj))
						return;
					else if (!is_dominated(temp_obj, y_obj))
						updatePA(archive);
				}
			}
		}
	}
}

void TIndividual::hybridNew(vector<vector<double>> &archive)
{
	if (rand_real(0, 1) < 0.5)
		swapNew(archive);
	else
		insertNew(archive);
}

void TIndividual::localIntensifization(vector<vector<double>> &archive)
{
	TIndividual sol0;
	sol0.x_var = x_var;
	sol0.x_speed = x_speed;
	sol0.y_obj = y_obj;
	bool costflag = 0;
	bool updateflag = 0;
	int cnt = 0, k = 0;
	int job = 0;
	int f;
	vector<pair<int, int>>jobPermu;

	for (int i = 0; i < factorysize; i++) {
		for (int j = 0; j < sol0.x_var[i].size(); j++)
			jobPermu.push_back(pair<int, int>{ i, sol0.x_var[i][j] });
	}
	while(k<jobsize && cnt<jobsize){
		job = jobPermu[k].second;
		f = jobPermu[k].first;
		if (sol0.x_var[f].size() <= 1) {
			k = (k + 1) % jobsize + 1;
			continue;
		}
		sol0.x_var[f].erase(find(sol0.x_var[f].begin(), sol0.x_var[f].end(),job));
		//sol0.factory_obj(f);
		if (costflag == 0) {
			sol0.randspeedup(f);
			sol0.speedup(f);
		}
		else {
			sol0.randspeeddown(f);
			sol0.speeddown(f);
		}
		vector<vector<int>>temp_speed(sol0.x_speed);
		for (int i = 0; i < factorysize; i++) {
			if (i != f) {
				for (int j = 0; j < sol0.x_var[i].size(); j++) {
					sol0.x_var[i].insert(sol0.x_var[i].begin() + j, job);
					//sol0.factory_obj(i);
					if (costflag == 0) {
						sol0.randspeedup(i);
						sol0.speedup(i);
					}
					else {
						sol0.randspeeddown(i);
						sol0.speeddown(i);
					}
					sol0.obj_eval();
					/*sol0.show_variable();
					sol0.show_objective();
					show_objective();
					cout << endl;*/
					if (is_dominated( sol0.y_obj, y_obj)) {
						x_speed = sol0.x_speed;
						x_var = sol0.x_var;
						y_obj = sol0.y_obj;
						updateflag = 1;
					}
					else if (!is_dominated(y_obj, sol0.y_obj)) {
						sol0.updatePA(archive);
					}
					sol0.x_speed = temp_speed;
					sol0.x_var[i].erase(sol0.x_var[i].begin() + j);
				}
			}
		}
		sol0.x_speed = x_speed;
		sol0.x_var = x_var;
		sol0.y_obj = y_obj;
		if (updateflag)
			cnt = 0;
		else
			cnt++;
		k = (k + 1) % jobsize + 1;
		updateflag = 0;
	}
}

void TIndividual::show_objective()
{
	for (int n = 0; n < numObjectives; n++)
		cout << y_obj[n] << " ";
	cout << endl;
}

void TIndividual::show_variable()
{
	for (vector<vector<int>>::iterator it = x_var.begin(); it != x_var.end(); it++) {
		for (vector<int>::iterator t = (*it).begin(); t != (*it).end(); t++) {
			cout << *t << " ";
		}
		cout << endl;
	}
}

void TIndividual::show_speed()
{
	for (vector<vector<int>>::iterator it = x_speed.begin(); it != x_speed.end(); it++) {
		for (vector<int>::iterator t = (*it).begin(); t != (*it).end(); t++) {
			cout << *t << " ";
		}
		cout << endl;
	}
}

int cmp1(const void* a, const void* b)
//目标函数f1的升序排序
{
	const TIndividual* e = (const TIndividual*)a;
	const TIndividual* f = (const TIndividual*)b;
	if (e->y_obj[0] == f->y_obj[0])
		return 0;
	else if (e->y_obj[0] < f->y_obj[0])
		return -1;
	else return 1;
}

int cmp2(const void* a, const void* b)
//目标函数f2的升序排序
{
	const TIndividual* e = (const TIndividual*)a;
	const TIndividual* f = (const TIndividual*)b;
	if (e->y_obj[1] == f->y_obj[1])
		return 0;
	else if (e->y_obj[1] < f->y_obj[1])
		return -1;
	else return 1;
}

int cmp_d(const void* a, const void* b)
//对拥挤距离降序排序
{
	const TIndividual* e = (const TIndividual*)a;
	const TIndividual* f = (const TIndividual*)b;
	if (e->crowding_distance == f->crowding_distance)
		return 0;
	else if (e->crowding_distance < f->crowding_distance)
		return 1;
	else
		return -1;
}
#endif
