#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
using namespace std;
vector<vector<int> > adj;
int Max(int a, int b, int c)
{
	int max;
	if (a >= b)
	{
		if (a >= c) {
			max = a;
		}
		else
			max = c;
	}
	else if (b >= c) { max = b; }
	else max = c;
	return max;
}
void Input_data(const string& filename) {
	std::fstream in(filename.c_str());
	cout << "reading file:" << filename << endl;
	string s;
	size_t n = 0, m = 0;
	string data1, data2;
	while (true)
	{
		std::getline(in, s);
		istringstream is(s);
		is >> data1 >> data2;
		int d1 = stoi(data1);
		int d2 = stoi(data2);
		n = Max(n, d2, d1);
		m += 1;
		if (in.eof()) { break; }
	}
	//this block will count the number of lines and calculate the maximun number appeared in the file, which are the parameters n, m(vertice, edge)
	in.seekg(0, ios::beg);
	n += 1;
	adj = vector<vector<int> >(n, vector<int>());
	for (size_t i = 0; i < m; i++) 
	{
		int x, y;
		std::getline(in, s);
		istringstream is(s);
		is >> data1 >> data2;
		x = stoi(data1);
		y = stoi(data2);
		adj[x].push_back(y);
	}
	in.close();
	//this block will assign data into the vertice template in terms of the adjancancy list
}

void Print_data() {
	int n = adj.size();
	for (int i = 0; i < n; i++) {
		int m = adj[i].size();
		cout << i << ": ";
		for (int j = 0; j < m; j++) { cout << adj[i][j] << " ,"; }
		cout << endl;

	}

}
int main() {
	Input_data("out.txt");
	Print_data();
}