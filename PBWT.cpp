#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <map>
#include <cassert>
#include <cstring>
#include <algorithm>

using namespace std;

int M, N, L, W, G;  // # of sequences, # of sites, length, width, gap size
double rho; // minimum allele frequency threshold

struct VCFReader {
	ifstream vcf;
	int G, M, p1 = 0; // p1 points to first site in gap
	vector<vector<int>> gap; // stores haplotype data in the gap
	vector<string> ID, fixedFields;

	VCFReader(string file, int _G, int _M, ofstream& out) {
		vcf = ifstream(file);
		G = max(_G, 1), M = _M; // gap size of 0 is equivalent to a gap size of 1 when handling the VCF file and gap
		gap = vector<vector<int>>(G, vector<int>(M));
		ID = vector<string>(M);
		fixedFields = vector<string>(G);
		preprocess(out);
		initGap();
	}

	VCFReader(string file, int _G, int _M) {
		vcf = ifstream(file);
		G = max(_G, 1), M = _M; // gap size of 0 is equivalent to a gap size of 1 when handling the VCF file and gap
		gap = vector<vector<int>>(G, vector<int>(M));
		ID = vector<string>(M);
		fixedFields = vector<string>(G);
		preprocess();
		initGap();
	}
	
	// gets haplotype IDs and moves input stream pointer to start of raw data
	void preprocess(ofstream& out) { 
		// skip meta-info lines and get header line
		string header;
		while (getline(vcf, header)) {
			out << header << '\n';
			if ((int)header.size() < 2 || header[0] != '#' || header[1] != '#') break;
		}

		// input sample IDs
		stringstream ss(header);
		for (int i = 0; i < 9; ++i) getline(ss, ID[0], '\t'); // skip fixed columns, assumes 9 columns (FORMAT column) 
		for (int i = 0; i < M / 2; ++i) {
			getline(ss, ID[2 * i], '\t');
			ID[2 * i + 1] = ID[2 * i] + "-1";
			ID[2 * i] += "-0";
		}
	}

	void preprocess() { 
		// skip meta-info lines and get header line
		string header;
		while (getline(vcf, header)) {
			if ((int)header.size() < 2 || header[0] != '#' || header[1] != '#') break;
		}

		// input sample IDs
		stringstream ss(header);
		for (int i = 0; i < 9; ++i) getline(ss, ID[0], '\t'); // skip fixed columns, assumes 9 columns (FORMAT column) 
		for (int i = 0; i < M / 2; ++i) {
			getline(ss, ID[2 * i], '\t');
			ID[2 * i + 1] = ID[2 * i] + "-1";
			ID[2 * i] += "-0";
		}
	}

	// initializes sliding window for the gap
	void initGap() {
		for (int i = 0; i < G; ++i) nextSite();
	}

	// reads the next site in the VCF file
	void nextSite() {
		fixedFields[p1].clear();
		char s[2 * M + 5000]; // assumes fixed fields take up less than 5000 characters
		vcf.getline(s, 2 * M + 5000);
		// skip fixed fields
		int offset = 0, tabs = 0; // offset = position in "s" of first sequence - points to the first character after 9 tabs
		while (tabs < 9) {
			if (s[offset] == '\t') ++tabs;
			if (tabs < 9) fixedFields[p1] += s[offset];
			++offset;
		}

		for (int i = 0; i < M; ++i) {
			assert(s[offset + (i / 2) * 4 + 1] == '|'); // sanity check
			gap[p1][i] = (s[offset + 2 * i] == '0' ? 0 : 1);
		}
		p1 = (p1 + 1) % G;
	}

	int getAllele(int g, int idx) {
		return gap[(p1 + g) % G][idx];
	}

	string getFixedField(int g) {
		return fixedFields[(p1 + g) % G];
	}

	void editGap(int g, int idx, int allele) {
		gap[(p1 + g) % G][idx] = allele;
	}

	void close() {vcf.close();}
};


struct biPBWT {
	int M, site;
	vector<int> pre, div, backwardPre, backwardDiv;
	vector<int> a, b, d, e;
	vector<int> block, rBlock;

	biPBWT() {}
	biPBWT(int _M) {
		this->M = _M;
		site = 0;
		pre = vector<int>(M);
		iota(pre.begin(), pre.end(), 0);
		div = vector<int>(M);
		backwardPre = vector<int>(M), backwardDiv = vector<int>(M);
		a = vector<int>(M), b = vector<int>(M), d = vector<int>(M), e = vector<int>(M);
		block = vector<int>(M), rBlock = vector<int>(M);
	}

	void nextSite(VCFReader& vcf, ifstream& backward) {
		// update forward PBWT
		int u = 0, v = 0, p = site + 1, q = site + 1;
		for (int i = 0; i < M; ++i) {
			int id = pre[i];
			if (div[i] > p) p = div[i];
			if (div[i] > q) q = div[i];
			if (vcf.getAllele(0, id) == 0) {
				a[u] = id;
				d[u] = p;
				++u;
				p = 0;
			}
			else {
				b[v] = id;
				e[v] = q;
				++v;
				q = 0;
			}
		}
		for (int i = 0; i < u; ++i) {
			pre[i] = a[i];
			div[i] = d[i];
		}
		for (int i = 0; i < v; ++i) {
			pre[u + i] = b[i];
			div[u + i] = e[i];
		}
		++site;

		// update reverse PBWT
		int rsite = (N - 1) - site - G; // index of the corresponding reverse site
		backward.seekg((long long)rsite * M * 8);

		// initialize backwardPre, backwardDiv, and rBlock
		int id = 0;
		for (int i = 0; i < M; ++i) {
			backward.read((char*)&backwardPre[i], sizeof backwardPre[i]);
			int rDiv; backward.read((char*)&rDiv, sizeof rDiv);
			backwardDiv[i] = rDiv;
			rDiv = (N - 1) - rDiv; // get forward index for position comparision

			if (rDiv < site + (G - 1) + L) ++id;
			rBlock[backwardPre[i]] = id;
		}

		// initialize block
		id = 0;
		for (int i = 0; i < M; ++i) {
			if (div[i] > site - L) ++id;
			block[pre[i]] = id;
		}

	}

	void countingSort(vector<vector<int>>& v, int idx) {
		vector<vector<vector<int>>> table(M + 1);
		for (int i = 0; i < M; ++i) {
			table[v[i][idx]].push_back(v[i]);
		}
		int p = 0;
		for (int i = 0; i <= M; ++i) {
			for (int j = 0; j < (int)table[i].size(); ++j) {
				v[p++] = table[i][j];
			}
		}
	}

	// gets blocks at the current site and processes them
	void getBlocks(VCFReader& orig, VCFReader& edit) {
		// Algorithm 2 - block matching
		vector<vector<int>> link(M, vector<int>(3)); // [sample ID, forward block ID, reverse block ID]
		for (int i = 0; i < M; ++i) {
			link[i][0] = i, link[i][1] = block[i], link[i][2] = rBlock[i];
		}
		// radix sort
		countingSort(link, 2);
		countingSort(link, 1);

		int start = 0;
		for (int i = 1; i < M; ++i) {
			if (link[i][1] != link[i - 1][1] || link[i][2] != link[i - 1][2]) {
				processBlock(link, start, i, orig, edit);
				start = i;	
			}
		}
		processBlock(link, start, M, orig, edit);
	}

	void processBlock(vector<vector<int>>& link, int start, int end, VCFReader& orig, VCFReader& edit) { // [start, end)
		int blockSize = end - start;
		if (blockSize < W) return; // width too small
		
		vector<int> zero(G), one(G);
		for (int j = start; j < end; ++j) {
			int id = link[j][0];
			for (int k = 0; k < G; ++k) {
				if (orig.getAllele(k, id) == 0) ++zero[k];
				else ++one[k];
			}
		}

		// error correct
		for (int k = 0; k < G; ++k) {
			if (min(zero[k], one[k]) <= blockSize * rho) {
				int correctAllele = zero[k] < one[k] ? 1 : 0;
				for (int j = start; j < end; ++j) {
					int id = link[j][0];
					edit.editGap(k, id, correctAllele);
				}
			}
		}
	}
};


int main(int argc, char* argv[]) {
	ios_base::sync_with_stdio(0); cin.tie(0);

	string writeTo = string(argv[2]);
	ifstream backward(writeTo + ".rpbwt"), meta(writeTo + ".meta"), gMap(argv[3]);
	ofstream out(writeTo + ".smooth.vcf");

	int checkpoint = atoi(argv[4]);
	L = atof(argv[5]), W = atoi(argv[6]), G = atoi(argv[7]), rho = atof(argv[8]);

	// retrieve M and N from meta file
	meta >> M >> N;

	VCFReader orig(string(argv[1]), G, M, out), edit(string(argv[1]), G, M); // orignal vcf and editted vcf

	biPBWT bipbwt(M);
	int site;
	for (site = 0; site + G < N; ++site) {
		out << edit.getFixedField(0);
		for (int i = 0; i < M / 2; ++i) {
			out << '\t' << edit.getAllele(0, 2 * i) << '|' << edit.getAllele(0, 2 * i + 1);
		}
		out << '\n';
		bipbwt.nextSite(orig, backward);

		orig.nextSite();
		edit.nextSite();

		bipbwt.getBlocks(orig, edit);

		if (site % checkpoint == 0) cout << "Checkpoint " << site << endl;
	}
	for (int j = 0; j < G; ++j) {
		out << edit.getFixedField(0);
		for (int i = 0; i < M / 2; ++i) {
			out << '\t' << edit.getAllele(j, 2 * i) << '|' << edit.getAllele(j, 2 * i + 1);
		}
		out << '\n';
		++site;
	}

	orig.close();
	edit.close();
	meta.close();
	gMap.close();
	out.close();

	return 0;
}
