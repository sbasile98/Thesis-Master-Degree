#include <bits/stdc++.h>
#include <random>
#include <getopt.h>

using namespace std;

/*-------------------------------------------------------------------------------------------------------------*/
float MUTATION_ODD;
float MUTATION_ODD_THRESHOLD = 25.0;
float MUTATION_ODD_MIN_LIMIT = 5.0;
float MUTATION_ODD_MAX_LIMIT = 100.0;

int THETA;
double LAMBDA;
double EPSILON;
double DELTA;
/*-------------------------------------------------------------------------------------------------------------*/

vector<string> split(string str, string sep) {
    char* cstr = const_cast<char*>(str.c_str());
    char* current;
    vector<string> arr;
    current = strtok(cstr,sep.c_str());
    while(current != NULL) {
        arr.push_back(current);
        current = strtok(NULL, sep.c_str());
    }
    return arr;
}

void readGraph(const string &filename, int *n, int *e, string &lw, string &uw, string &seed, vector<int> &w, vector<vector<int> > &adj) {
    ifstream in(filename);
    string s;
    for (int i = 0; i < 3; i++) getline(in, s);

    in >> s;
    in >> *n;
    in >> s;
    in >> *e;
    in >> s;
    in >> s;

    lw = split(s, "-")[0];
    uw = split(s, "-")[1];

    in >> s;
    in >> seed;

    for (int i = 0; i < 2; i++) getline(in, s);

    for (int i = 0; i < *n; i++) {
        int a;
        in >> a;
        in >> a;
        w.push_back(a);
    }

    getline(in, s);
    getline(in, s);

    adj.resize(*n);
    for (int i = 0; i < *n; i++) adj[i].resize(*n);

    for (int i = 0; i < *n; i++) {
        for (int j = 0; j <= i; j++) {
            int a;
            in >> a;

            adj[i][j] = a;
            adj[j][i] = a;
        }
    }
}

void printInput(int n, int e, vector<int> &w, vector<vector<int> > &adj) {
    cout << n << " " << e << endl;

    for (int i = 0; i < n; i++) {
        cout << w[i] << " ";
    }
    cout << endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            cout << adj[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int operator==(vector<int> &x, vector<int> &y) {
    for (int i = 0; i < x.size() - 2; i++) {
        if (x[i] != y[i]) return 0;
    }
    return 1;
}

int fitness(vector<int> &x, vector<int> &w, int n) {
    int y = 0;

    for (int i = 0; i < n; i++) {
        y += x[i] * w[i];
    }

    return y;
}

int thetaf(vector<int> &w, int n) {
    int penalty = 0;

    for (int i = 0; i < n; i++) penalty += w[i];

    return penalty;
}

int DFS(vector<vector<int> > &adj, vector<int> &x, int u, int *color, int *p, int n) {
    color[u] = 1;
    int res = 0;

    for (int j = 0; j < n; j++) {
        if (adj[u][j] && p[u] != j && !x[u] && !x[j]) {
            if (color[j] == 1) {
                return 1;
            }
            if (color[j] == 0) {
                p[j] = u;
                res += DFS(adj, x, j, color, p, n);
                if (res > 0) return 1;
            }
        }
    }

    color[u] = 2;
    return res;
}

int DFS(vector<vector<int> > &adj, vector<int> &x, int u, int *color, int *p, int n, vector<int> &conn) {
    color[u] = 1;
    int res = 0;

    conn[u] = 1;

    for (int j = 0; j < n; j++) {
        if (adj[u][j] && p[u] != j && !x[u] && !x[j]) {
            if (color[j] == 1) {
                return 1;
            }
            if (color[j] == 0) {
                p[j] = u;
                res += DFS(adj, x, j, color, p, n, conn);
                if (res > 0) return 1;
            }
        }
    }

    color[u] = 2;
    return res;
}

int hasCycle(vector<int> &x, vector<vector<int> > &adj, int n) {
    int y = 0;

    int color[n];
    for (int i = 0; i < n; i++) color[i] = 0;

    int p[n];
    for (int i = 0; i < n; i++) p[i] = -1;

    for (int i = 0; i < n; i++) {
        if (color[i] == 0) y += DFS(adj, x, i, color, p, n);
    }

    return y;
}

int hasCycle(vector<int> &x, vector<vector<int> > &adj, int n, int u) {
    int y = 0;

    int color[n];
    for (int i = 0; i < n; i++) color[i] = 0;

    int p[n];
    for (int i = 0; i < n; i++) p[i] = -1;

    y += DFS(adj, x, u, color, p, n);

    return y;
}

void connectedComponents(vector<int> &x, vector<vector<int> > &adj, int n, int u, vector<int> &conn) {
    for (int i = 0; i < n; i++) conn[i] = 0;

    int color[n];
    for (int i = 0; i < n; i++) color[i] = 0;

    int p[n];
    for (int i = 0; i < n; i++) p[i] = -1;

    DFS(adj, x, u, color, p, n, conn);
    conn[0] = 0;
}

int heavierVertex(vector<int> &x, vector<int> &w, int n) {
    int max = numeric_limits<int>::min();
    int index = -1;

    for (int i = 0; i < n; i++) {
        if (x[i] && w[i] > max) {
            max = w[i];
            index = i;
        }
    }

    return index;
}

int vertexDegree(int u, vector<int> &x, vector<vector<int> > &adj, int n) {
    int degree = 0;

    for (int i = 0; i < n; i++) {
        if (!x[i]) degree += adj[u][i];
    }

    return degree;
}

int minWeightDegreeRatio(int u, vector<int> &x, vector<int> &w, vector<vector<int> > &adj, int n, vector<int> &conn) {
    float min = numeric_limits<float>::max();
    int index = -1;

    for (int i = 0; i < n; i++) {
        if (!x[i] && conn[i]) {
            float ratio = w[i] / ceil(vertexDegree(i, x, adj, n));
            if (ratio < min) {
                min = ratio;
                index = i;
            }
        }
    }

    return index;
}

int maxWeightDegreeRatio(int u, vector<int> &x, vector<int> &w, vector<vector<int> > &adj, int n) {
    float max = numeric_limits<float>::min();
    int index = -1;

    for (int i = 0; i < n; i++) {
        if (x[i]) {
            x[i] = 0;
            float ratio = w[i] / ceil(vertexDegree(i, x, adj, n));
            if (ratio > max) {
                max = ratio;
                index = i;
            }
            x[i] = 1;
        }
    }

    return index;
}

int solutionCardinality(vector<int> &x, int n) {
    int count = 0;

    for (int i = 0; i < n; i++) if(x[i]) count++;

    return count;
}

void mutation(vector<int> &x, vector<vector<int> > &adj, int n, int mutationType) {
    random_device rd;

    if (mutationType == 0) {
        for (int i = 0; i < n; i++) {
            int odd = rd() % 100;

            if (odd < MUTATION_ODD && MUTATION_ODD < MUTATION_ODD_THRESHOLD) x[i] = 1 - x[i];
            if (odd < MUTATION_ODD && MUTATION_ODD >= MUTATION_ODD_THRESHOLD && !x[i]) x[i] = 1 - x[i];
        }
    }
    else if (mutationType == 1) {
        for (int i = 0; i < n; i++) {
            int odd = (((rd() % 100) + (vertexDegree(i, x, adj, n) / n * 100)) / 2);

            if (odd < MUTATION_ODD && MUTATION_ODD < MUTATION_ODD_THRESHOLD && !x[i]) x[i] = 1 - x[i];
            if (odd < MUTATION_ODD && MUTATION_ODD >= MUTATION_ODD_THRESHOLD && !x[i]) x[i] = 1 - x[i];
        }
    }
}

vector<vector<int> > crossover(vector<int> &p_1, vector<int> &p_2, int n, int type) {
    vector<vector<int> > s;
    random_device rd;

    s.resize(2);
    for (int i = 0; i < 2; i++) s[i].resize(n + 2);

//    1-point crossover
    if (type == 1) {
        int split = rd() % (n - 2);

        for (int i = 0; i < n; i++) {
            if (i <= split) {
                s[0][i] = p_1[i];
                s[1][i] = p_2[i];
            }
            else {
                s[0][i] = p_2[i];
                s[1][i] = p_1[i];
            }
        }
    }

//    2-point crossover
    else if (type == 2) {
        int split_1 = rd() % (n - 3);
        int split_2 = (split_1 + 1) + (rd() % (n - split_1 - 2));

        for (int i = 0; i < n; i++) {
            if (i <= split_1) {
                s[0][i] = p_1[i];
                s[1][i] = p_2[i];
            }
            else if (i > split_1 && i <= split_2){
                s[0][i] = p_2[i];
                s[1][i] = p_1[i];
            }
            else {
                s[0][i] = p_1[i];
                s[1][i] = p_2[i];
            }
        }
    }

//    3-point crossover
    else if (type == 3) {
        int split_1 = rd() % (n - 4);
        int split_2 = (split_1 + 1) + (rd() % (n - split_1 - 3));
        int split_3 = (split_2 + 1) + (rd() % (n - split_2 - 2));

        for (int i = 0; i < n; i++) {
            if (i <= split_1) {
                s[0][i] = p_1[i];
                s[1][i] = p_2[i];
            }
            else if (i > split_1 && i <= split_2){
                s[0][i] = p_2[i];
                s[1][i] = p_1[i];
            }
            else if (i > split_2 && i <= split_3){
                s[0][i] = p_1[i];
                s[1][i] = p_2[i];
            }
            else {
                s[0][i] = p_2[i];
                s[1][i] = p_1[i];
            }
        }
    }

//    Uniform crossover
    else if (type == 4) {
        for (int i = 0; i < n; i++) {
            s[0][i] = rd() % 2;
            s[1][i] = 1 - s[0][i];
        }
    }

    return s;
}

void createPopulation(int n, int m, vector<vector<int> > &pop, vector<vector<int> > &adj, vector<int> &w) {
    random_device rd;
    int max_w = *max_element(w.begin(), w.end());

    pop.resize(m);
    for (int i = 0; i < m; i++) pop[i].resize(n + 2);

//    Effettuo la costruzione della popolazione iniziale con un approccio stocastico misto ad uno deterministico facendo in modo che un nodo venga inserito in soluzione con una probabilità direttamente proporzionale al suo peso
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            pop[i][j] = (((rd() % 100) + ((double(w[j]) / double(max_w)) * 100)) / 2) <= 50 ? 1 : 0;
        }
    }

    for (int i = 0; i < m; i++) {
        pop[i][n] = hasCycle(pop[i], adj, n) ? (THETA + THETA/ceil(fitness(pop[i], w, n))) : fitness(pop[i], w, n);
        pop[i][n + 1] = hasCycle(pop[i], adj, n) ? 0 : 1;
    }
}

void printPopulation(vector<vector<int> > &pop, int n, int m) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n + 2; j++) {
            cout << pop[i][j] << " ";
            if (j == n - 1) cout << " | ";
        }
        cout << endl;
    }
    cout << endl << endl;
}

bool sortcol(const vector<int> &v1, const vector<int> &v2) {
    return v1[v1.size() - 2] < v2[v1.size() - 2];
}

vector<vector<int> > createChildren(vector<vector<int> > &pop, int n, int m, vector<vector<int> > &adj, vector<int> &w, int crossoverType, int mutationType) {
    random_device rd;
    vector<vector<int> > children;

    for (int i = 0; i < m/2; i++) {
        int a = rd() % m;
        int b = rd() % m;

        while (pop[a] == pop[b]) b = rd() % m;

        vector<int> p_1 = pop[a];
        vector<int> p_2 = pop[b];

        vector<vector<int> >tmp_child = crossover(p_1, p_2, n, crossoverType);
        for (int j = 0; j < 2; j++) mutation(tmp_child[j], adj, n, mutationType);

        for (int j = 0; j < 2; j++) {
            tmp_child[j][n] = hasCycle(tmp_child[j], adj, n) ? (THETA + THETA/ceil(fitness(tmp_child[j], w, n))) : fitness(tmp_child[j], w, n);
            tmp_child[j][n + 1] = hasCycle(tmp_child[j], adj, n) ? 0 : 1;

            if (fitness(tmp_child[j], w, n) == 0) tmp_child[j][n] = THETA * 2;
        }

        for (int j = 0; j < 2; j++) children.push_back(tmp_child[j]);
    }

    return children;
}

void printChildren(vector<vector<int> > &children, int &n, int &m) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n + 2; j++) {
            cout << children[i][j] << " ";
            if (j == n - 1) cout << " | ";
        }
        cout << endl;
    }
    cout << endl << endl;
}

void nextGeneration(vector<vector<int> > &pop, vector<vector<int> > &children, vector<vector<int> > &adj, int n, int m, int joinType) {
    if (joinType == 0) {
        for (int i = 0; i < m; i++) pop.push_back(children[i]);
        sort(pop.begin(), pop.end(), sortcol);

        for (int i = 0; i < m; i++) {
            pop.pop_back();
        }
    }
    else if (joinType == 1) {
        sort(pop.begin(), pop.end(), sortcol);
        sort(children.begin(), children.end(), sortcol);

        for (int i = 0; i < m/2; i++) {
            pop.pop_back();
        }
        for (int i = 0; i < m/2; i++) {
            pop.push_back(children[i]);
        }
    }
    if (joinType == 2) {
        for (int i = 0; i < m; i++) pop.push_back(children[i]);
        sort(pop.begin(), pop.end(), sortcol);

        int i = 1;
        int rem = 0;

        while (i < pop.size() && rem < m) {
            if (pop[i] == pop[i - 1]) {
                pop.erase(pop.begin() + i);
                rem++;
            }
            i++;
        }

        for (int i = rem; i < m; i++) pop.pop_back();
    }
}

/*-------------------------------------------------------------------------------------------------------------*/
double meanFitness(vector<vector<int> > &pop, int n, vector<int> &w) {
    double sumFitness = 0;

    for (int i = 0; i < n; i++) sumFitness += fitness(pop[i], w, n);

    return sumFitness / n;
}

double meanBufferFitness(deque<double> &fitnessBuffer, int v) {
    double sumFitness = 0;

    for (int i = 0; i < v; i++) sumFitness += fitnessBuffer[i];

    return sumFitness / v;
}

double correlation(deque<double> &fitnessBuffer, int v) {
    double meanFitness = meanBufferFitness(fitnessBuffer, v);

    double corr = 0;
    for (int i = 0; i < v; i++) corr += (fitnessBuffer[i] - meanFitness) * (i - v/2 - 0.5);

    return corr;
}

void computeMutationOdd(deque<double> &fitnessBuffer, int v) {
	double dev_stand = (pow(v, 2) - 1)  / 12;
    double corr = correlation(fitnessBuffer, v);
	double alpha = corr / dev_stand;

    if (alpha > EPSILON) MUTATION_ODD -= DELTA;
    else if (abs(alpha) < EPSILON) MUTATION_ODD += DELTA * LAMBDA; 
    else if (alpha < -EPSILON) MUTATION_ODD += DELTA;
	
	if (MUTATION_ODD < MUTATION_ODD_MIN_LIMIT) MUTATION_ODD = MUTATION_ODD_MIN_LIMIT;
	if (MUTATION_ODD > MUTATION_ODD_MAX_LIMIT) MUTATION_ODD = MUTATION_ODD_MAX_LIMIT;	
}

void mutation(vector<int> &x, vector<vector<int> > &adj, int n, deque<double> &fitnessBuffer, int v) {
    random_device rd;
	
	double dev_stand = (pow(v, 2) - 1)  / 12;
    double corr = correlation(fitnessBuffer, v);
	double alpha = corr / dev_stand;
	
	for (int i = 0; i < n; i++) {
		int odd = (((rd() % 100) + (vertexDegree(i, x, adj, n) / n * 100)) / 2);

		if (odd < MUTATION_ODD && abs(alpha) < EPSILON && !x[i]) x[i] = 1 - x[i];
		else if (odd < MUTATION_ODD) x[i] = 1 - x[i];
	}
}

vector<vector<int> > createChildren(vector<vector<int> > &pop, int n, int m, vector<vector<int> > &adj, vector<int> &w, int crossoverType, deque<double> &fitnessBuffer, int v) {
    random_device rd;
    vector<vector<int> > children;

    for (int i = 0; i < m/2; i++) {
        int a = rd() % m;
        int b = rd() % m;

        while (pop[a] == pop[b]) b = rd() % m;

        vector<int> p_1 = pop[a];
        vector<int> p_2 = pop[b];

        vector<vector<int> >tmp_child = crossover(p_1, p_2, n, crossoverType);
        for (int j = 0; j < 2; j++) mutation(tmp_child[j], adj, n, fitnessBuffer, v);

        for (int j = 0; j < 2; j++) {
            tmp_child[j][n] = hasCycle(tmp_child[j], adj, n) ? (THETA + THETA/ceil(fitness(tmp_child[j], w, n))) : fitness(tmp_child[j], w, n);
            tmp_child[j][n + 1] = hasCycle(tmp_child[j], adj, n) ? 0 : 1;

            if (fitness(tmp_child[j], w, n) == 0) tmp_child[j][n] = THETA * 2;
        }

        for (int j = 0; j < 2; j++) children.push_back(tmp_child[j]);
    }

    return children;
}
/*-------------------------------------------------------------------------------------------------------------*/

void localSearch(vector<vector<int> > &pop, vector<vector<int> > &adj, vector<int> &w, int n, int inf, int sup, int &min) {
    for (int i = inf; i < sup; i++) {
        vector<int> copy = pop[i];
        int it = solutionCardinality(pop[i], n);

        for (int j = 0; j < it; j++) {
//            Seleziono il nodo con rapporto peso/grado massimo, lo rimuovo dalla soluzione e dalla copia di quest'ultima
            int u = maxWeightDegreeRatio(u, copy, w, adj, n);
            pop[i][u] = 0;
            copy[u] = 0;

            int v;
            bool find = false;

//            Se è presente un ciclo nel sottografo a partire dal nodo u che ho rimosso, cerco e provo a rimuovere un nodo v con rapporto peso/grado minimo
            if (hasCycle(pop[i], adj, n, u)) {
                vector<int> conn;
                conn.resize(n);
                connectedComponents(pop[i], adj, n, u, conn);

//                Per ogni componenente connessa di u eseguo la ricerca di un nodo v di rapporto peso/grado minimo la cui rimozione rende aciclico il sottografo
                for (int k = 0; k < n; k++) {
                    v = minWeightDegreeRatio(u, pop[i], w, adj, n, conn);
                    conn[v] = 0;
                    pop[i][v] = 1;

                    if (!hasCycle(pop[i], adj, n, u)) {
                        find = true;
                        break;
                    }
                    else {
                        pop[i][v] = 0;
                    }
                }

//                Se ho trovato un nodo il cui inserimento in soluzione rende il grafo aciclico lo rimuovo dal sottografo altrimenti reinserisco in soluzione il nodo u
                if (find) {
                    pop[i][n] = pop[i][n] - w[u] + w[v];
                    find = false;
                }
                else {
                    pop[i][u] = 1;
                }
            }
            else {
                pop[i][n] -= w[u];
            }
        }
    }

    sort(pop.begin(), pop.end(), sortcol);
    if (min > pop[0][n]) min = pop[0][n];
}

int main(int argc, char *argv[]) {
	
	int opt;
	string filename;
	
	while ((opt = getopt(argc, argv, "f:e:d:l:m:")) != -1) {
		switch(opt) {
			case 'f':
				filename = optarg;
				break;
			
			case 'e':
				EPSILON = atof(optarg);
				break;
				
			case 'd':
				DELTA = atof(optarg);
				break;
				
			case 'l':
				LAMBDA = atof(optarg);
				break;
				
			case 'm':
				MUTATION_ODD = atof(optarg);
				break;
		}
	}
	
    int n, e, m = 100, crossoverType = 1, mutationType = 1, joinType = 2;
    string lw, uw, seed;
    int min = numeric_limits<int>::max(), count = 0, generationNumber = 1000;
    vector<int> w;
    vector<vector<int> > pop, adj, children;
	
	deque<double> fitnessBuffer;
	int v = 20;

    readGraph(filename, &n, &e, lw, uw, seed, w, adj);
    THETA = thetaf(w, n);

//    Creo la popolazione iniziale
    createPopulation(n, m, pop, adj, w);
	
    for (int i = 0; i < generationNumber; i++) {
//        Creo la popolazione dei figli applicando gli operatori di crossover e mutation
        if (i < v) children = createChildren(pop, n, m, adj, w, crossoverType, mutationType);
		else children = createChildren(pop, n, m, adj, w, crossoverType, fitnessBuffer, v);

//        Costruisco la popolazione della generazione successiva
        nextGeneration(pop, children, adj, n, m, joinType);

//        Local Search: cerco delle soluzioni migliori sulla totalità della popolazione
        localSearch(pop, adj, w, n, 0, m, min);
		
		if (min > pop[0][n] && pop[0][n] != 0) min = pop[0][n];
		
/*----------------------------------------------------------------------------------------------------------------------*/			
		if (i < v) {
            fitnessBuffer.push_back(meanFitness(pop, n, w));
        }
        else {
            fitnessBuffer.pop_front();
            fitnessBuffer.push_back(meanFitness(pop, n, w));
            computeMutationOdd(fitnessBuffer, v);
        }
/*----------------------------------------------------------------------------------------------------------------------*/
    }

//    Rimozione dei nodi ridondanti dalle soluzioni della popolazione finale
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            if (pop[i][j]) {
                pop[i][j] = 0;

                if (hasCycle(pop[i], adj, n, j)) pop[i][j] = 1;
                else {
                    pop[i][n] -= w[j];
                }
            }
        }
    }
    sort(pop.begin(), pop.end(), sortcol);
    if (min > pop[0][n]) min = pop[0][n];

//    Local Search: cerco delle soluzioni migliori sulla totalità della popolazione
    localSearch(pop, adj, w, n, 0, m, min);

//    Costruzione report dei risultati
    string fvs = "";

    for (int i = 0; i < n; i++) {
        string str = to_string(pop[0][i]);
        fvs += str;
    }

    ofstream out("results.csv", std::ofstream::app);
    out << filename << "\t" << n << "\t" << e << "\t" << lw << "\t" << uw << "\t" << seed << "\t" << min << "\t" << fvs << endl;

    cout << filename << "\t" << min << "\t" << (pop[0][n + 1] ? ("Feasible") : ("Unfeasible")) << "\t" << fvs << endl;

    return 0;
}
