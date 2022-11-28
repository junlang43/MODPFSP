#include "insga2.h"

int main()
{
	srand((unsigned int)time(NULL));
	for (size_t k = 1; k <= 1; k++) {
		ofstream outfile;
		string filename = "Perato20_16_4_" + to_string(k) + ".txt";
		outfile.open(filename, ostream::app);
		ReadFile(k);
		vector<vector<double>>archive;
		for (int i = 0; i < 10; i++) {
			TINSGA* pop = new TINSGA;
			pop->run();
			bool flag;
			for (auto it1 = pop->PA.begin(); it1 != pop->PA.end(); it1++) {
				flag = true;
				for (vector<vector<double>>::iterator it2 = archive.begin(); it2 != archive.end();) {
					if (is_dominated(*it1, *it2))
						it2 = archive.erase(it2);
					else {
						if (((*it1)[0] == (*it2)[0] && (*it1)[1] == (*it2)[1]) || is_dominated(*it2, *it1)) {
							flag = false;
							break;
						}
						it2++;
					}
				}
				if (flag) {
					archive.push_back(*it1);
				}
			}
		}
		for (auto it = archive.begin(); it != archive.end(); it++)
			outfile << (*it)[0] << " " << (*it)[1] << endl;
	}
	
	
	return 0;
}