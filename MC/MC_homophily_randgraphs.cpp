#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include "math.h"

using namespace std;

void init_ran(gsl_rng * &r, unsigned long s){
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, s);
}


typedef struct{
    vector <long> edges_in; // edges that contain the node
    vector <int> opinion;
    vector <long> neighs;
}Tnode;




typedef struct{
    vector <long> nodes_in; // nodes inside the edge. nodes[i], with i={0, 1}.
    int dot_prod;
    int dot_prod_new;
}Tedge;




void init_graph_from_file(Tnode *&nodes, Tedge *&edges, long &N, long &M, int g, char *graphname){
    ifstream fin(graphname);
    fin >> N;
    fin >> M;
    nodes = new Tnode[N];
    edges = new Tedge[M];
    long i, j;
    for (long e = 0; e < M; e++){
        fin >> i;
        fin >> j;
        edges[e].nodes_in.push_back(i - 1);
        edges[e].nodes_in.push_back(j - 1);
        

        nodes[i - 1].edges_in.push_back(e);
        nodes[j - 1].edges_in.push_back(e);
        nodes[i - 1].neighs.push_back(j - 1);
        nodes[j - 1].neighs.push_back(i - 1);
    }

    for (long nod = 0; nod < N; nod++){
        nodes[nod].opinion = vector <int> (g);
    }
    fin.close();
}


void rand_init(Tnode *nodes, long N, int g, gsl_rng * r, double p){
    for (long i = 0; i < N; i++){
        for (long l = 0; l < g; l++){
            nodes[i].opinion[l] = 1 - 2 * (gsl_rng_uniform(r) < p);
        }
    }
}

void ord_init(Tnode *nodes, long N, int g){
    for (long i = 0; i < N; i++){
        for (long l = 0; l < g; l++){
            nodes[i].opinion[l] = 1;
        }
    }
}


void init_dot_prod(Tnode *nodes, Tedge *edges, long M, int g){
    for (long e = 0; e < M; e++){
        edges[e].dot_prod = 0;
        for (long l = 0; l < g; l++){
            edges[e].dot_prod += nodes[edges[e].nodes_in[0]].opinion[l] * 
                                nodes[edges[e].nodes_in[1]].opinion[l];
        }
    }
}


double deltaE(Tnode *nodes, Tedge *edges, int pos_node, int component, double alpha, int g){
    int si_l, sj_l;
    double dE = 0;
    for (int j = 0; j < nodes[pos_node].neighs.size(); j++){
        si_l = nodes[pos_node].opinion[component];
        sj_l = nodes[nodes[pos_node].neighs[j]].opinion[component];
        edges[nodes[pos_node].edges_in[j]].dot_prod_new = edges[nodes[pos_node].edges_in[j]].dot_prod - 2 * si_l * sj_l;
        dE += (2 * alpha - 1) * 2 * si_l * sj_l - 
               fabs(edges[nodes[pos_node].edges_in[j]].dot_prod_new) + fabs(edges[nodes[pos_node].edges_in[j]].dot_prod);
    }
    return dE * 0.5 / g;
}


void propose_flip(Tnode *nodes, Tedge *edges, long N, double alpha, int g, double temp, gsl_rng * r){
    long pos_node = gsl_rng_uniform_int(r, N);
    int component = gsl_rng_uniform_int(r, g);
    double dE = deltaE(nodes, edges, pos_node, component, alpha, g);
    if (dE <= 0){
        nodes[pos_node].opinion[component] *= -1;
        for (int j = 0; j < nodes[pos_node].neighs.size(); j++){
            edges[nodes[pos_node].edges_in[j]].dot_prod = edges[nodes[pos_node].edges_in[j]].dot_prod_new;
        }
    }else{
        double p = exp(-dE / temp);
        if (gsl_rng_uniform(r) < p){
            nodes[pos_node].opinion[component] *= -1;
            for (int j = 0; j < nodes[pos_node].neighs.size(); j++){
                edges[nodes[pos_node].edges_in[j]].dot_prod = edges[nodes[pos_node].edges_in[j]].dot_prod_new;
            }
        }
    }
}


double magnetization(Tnode *nodes, long N, int g){
    vector <double> m(g);
    for (long i = 0; i < N; i++){
        for (int l = 0; l < g; l++){
            m[l] += nodes[i].opinion[l];
        }
    }
    double mag = 0;
    for (int l = 0; l < g; l++){
        m[l] /= N;
        mag += m[l] * m[l];
    }
    return sqrt(mag / g);
}


double energy(Tedge *edges, long M, double alpha, int g){
    double E = 0;
    for (long e = 0; e < M; e++){
        E -= edges[e].dot_prod * (2 * alpha - 1) + fabs(edges[e].dot_prod);
    }
    return 0.5 * E / M / g;
}

void MonteCarloHist(Tnode *nodes, Tedge *edges, long N, long M, double alpha, int g, double temp, 
                    long tlimit, long save_every, 
                    vector <double> &mag, vector <double> &ener, string string_init, double p,
                    gsl_rng * r){
    
    if (string_init == "rand"){
        rand_init(nodes, N, g, r, p);
    }else if (string_init == "ord"){
        ord_init(nodes, N, g);
    }else{
        cout << "Error: init string not recognized" << endl;
        exit(1);
    }

    init_dot_prod(nodes, edges, M, g);

    double m, E;
    m = magnetization(nodes, N, g);
    E = energy(edges, M, alpha, g);
    long counter = 0;
    mag[counter] += m;
    ener[counter] += E;
    counter++;
    cout << 0 << "\t" << m << "\t" << E << endl;
    for (long t = 0; t < tlimit; t++){
        for (long i = 0; i < N * g; i++){
            propose_flip(nodes, edges, N, alpha, g, temp, r);
        }
        if ((t + 1) % save_every == 0){
            m = magnetization(nodes, N, g);
            E = energy(edges, M, alpha, g);
            mag[counter] += m;
            ener[counter] += E;
            counter++;
            cout << t + 1 << "\t" << m << "\t" << E << endl;
        }
    }

}


void MonteCarlo_OneGraph(double alpha, int g, double temp, long tlimit, long save_every, 
                         long seed0, long nsamples, char *graphname, vector <double> &mag, 
                         vector <double> &ener, string string_init, double p){
    Tnode *nodes;
    Tedge *edges;
    long N, M;
    init_graph_from_file(nodes, edges, N, M, g, graphname);
    gsl_rng * r;
    init_ran(r, seed0); 
    for (long i = 0; i < nsamples; i++){
        cout << endl;
        cout << "# MC seed: " << seed0 + i << endl;
        MonteCarloHist(nodes, edges, N, M, alpha, g, temp, tlimit, save_every, 
                       mag, ener, string_init, p, r);
    }
    delete [] nodes;
    delete [] edges;
    gsl_rng_free(r);
}

void print_results(vector <double> mag, vector <double> ener, long ngraphs, long nsamples, 
                   char *filemag, char *fileener, long save_every){
    ofstream fmag(filemag);
    ofstream fener(fileener);
    fmag << "# ngraphs=" << ngraphs << "  nhist=" << nsamples << endl;
    fener << "# ngraphs=" << ngraphs << "  nhist=" << nsamples << endl;
    for (long i = 0; i < mag.size(); i++){
        fmag << i * save_every << "\t" << mag[i] / ngraphs / nsamples << endl;
        fener << i * save_every << "\t" << ener[i] / ngraphs / nsamples << endl;
    }
    fmag.close();
    fener.close();
}

void MonteCarlo_All(double alpha, int g, double temp, long tlimit, long save_every, long print_every,
                    long seedgraph0, long ngraphs, long seedhist0, long nsamples, 
                    char *graph_str0, string string_init, char *filemag, char *fileener, 
                    double p=0.5){
    long length = tlimit / save_every + 1;
    vector <double> mag(length);
    vector <double> ener(length);
    char graphname[200];
    long seedgraph;
    for (long seedgraph = seedgraph0; seedgraph < seedgraph0 + ngraphs; seedgraph++){
        sprintf(graphname, "%s_seed_%li.txt", graph_str0, seedgraph);
        cout << endl;
        cout << "# Graph: " << graphname << endl;
        cout << endl;
        MonteCarlo_OneGraph(alpha, g, temp, tlimit, save_every, seedhist0, nsamples, graphname, mag, ener, string_init, p);
        if ((seedgraph - seedgraph0) % print_every == 0){
            print_results(mag, ener, seedgraph - seedgraph0 + 1, nsamples, filemag, fileener, save_every);
        }
    }
}


int main(int argc, char *argv[]) {
    double alpha = atof(argv[1]);
    int g = atoi(argv[2]);
    double temp = atof(argv[3]);
    long tlimit = atol(argv[4]);
    long save_every = atol(argv[5]);
    long print_every = atol(argv[6]);
    long seedgraph0 = atol(argv[7]);
    long ngraphs = atol(argv[8]);
    long seedhist0 = atol(argv[9]);
    long nsamples = atol(argv[10]);
    char *graph_str0 = argv[11];
    string string_init = argv[12];
    char * graph_id = argv[13];

    char filemag[200];
    char fileener[200];

    sprintf(filemag, "MC_%s_mag_alpha_%.3lf_g_%d_temp_%.3lf_tl_%li_tsampl_%li_sgraph0_%li_shist0_%li_nhist_%li.txt",
            graph_id, alpha, g, temp, tlimit, save_every, seedgraph0, seedhist0, nsamples);
    sprintf(fileener, "MC_%s_ener_alpha_%.3lf_g_%d_temp_%.3lf_tl_%li_tsampl_%li_sgraph0_%li_shist0_%li_nhist_%li.txt",
            graph_id, alpha, g, temp, tlimit, save_every, seedgraph0, seedhist0, nsamples);

    MonteCarlo_All(alpha, g, temp, tlimit, save_every, print_every, seedgraph0, ngraphs, seedhist0, nsamples, graph_str0, string_init, filemag, fileener);    
    return 0;
}