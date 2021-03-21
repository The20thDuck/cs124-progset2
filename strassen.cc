#include <cstdio>
#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <cmath>
#include <cstdlib>
#include <vector> 
#include <time.h>
#define MAXN 2000


int stras[MAXN];
int reg[MAXN];

void task1() {
    stras[1] = 1;
    reg[1] = 1;
    for (int i = 2; i < MAXN; i++) {
        int c = (i+1)/2;
        int f = i/2;
        stras[i] = 7*stras[(i+1)/2] + 18*c*c;
        reg[i] = i*i*(2*i-1);
        printf("%d: %d %d\n", i, stras[i], reg[i]);
        if (stras[i] < reg[i]){
            break;
        }
    }
}




int numBits(int n) { 
    int num = 0;
    while (n) {
        n >>= 1;
        num++;
    }  
    return num;
}

#define NUM_INTERMED 10
#define NUM_P 7
int num_layers;
std::vector<std::vector<int>> med[NUM_INTERMED]; // med[i][j] stores the ith intermediate array with numBits(d)-1 = j 
std::vector<std::vector<int>> p[NUM_P]; // p[i][j] stores the ith p array with numBits(d)-1 = j 

// Adds a*b and stores the result in c. 
void add(int d, std::vector<int> &a, std::vector<int> &b, std::vector<int> &c) {
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            c[d*i + j] = a[d*i + j] + b[d*i + j];
        }
    }
}

// Subtracts a*b and stores the result in c. 
void sub(int d, std::vector<int> &a, std::vector<int> &b, std::vector<int> &c) {
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            c[d*i + j] = a[d*i + j] - b[d*i + j];
        }
    }
}

void printVec(int d, std::vector<int> &v) {
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            printf("%d ", v[d*i + j]);
        }
        printf("\n");
    }
}


// Multiply a*b and store the result in c
void strassen(int d, std::vector<int> &a, std::vector<int> &b, std::vector<int> &c) {
    // printf("stras %d\n", d);
    if (d == 1){
        c[0] = a[0]*b[0];
        return; 
    }
    else {
        int d2 = (d+1)/2;
        int layer = numBits(d2)-1;
        // printf("d2 %d, layer %d\n", d2, layer);
        for (int i = 0; i < d2; i++) {
            for (int j = 0; j < d2; j++) {

                med[0][layer][d2*i + j] = a[d*i + j]; // A
                if (j+d2 < d) 
                    med[1][layer][d2*i + j] = a[d*i + j+d2]; // B
                if (i+d2 < d)
                    med[2][layer][d2*i + j] = a[d*(i+d2) + j]; // C
                if (i+d2 < d && j+d2 < d)
                    med[3][layer][d2*i + j] = a[d*(i+d2) + j+d2]; // D
                med[4][layer][d2*i + j] = b[d*i + j]; // E
                if (j+d2 < d)
                    med[5][layer][d2*i + j] = b[d*i + j+d2]; // F
                if (i+d2 < d)
                    med[6][layer][d2*i + j] = b[d*(i+d2) + j]; // G
                if (i+d2 < d && j+d2 < d)
                    med[7][layer][d2*i + j] = b[d*(i+d2) + j+d2]; // H
            }
        }

        sub(d2, med[5][layer], med[7][layer], med[8][layer]); // F-H
        strassen(d2, med[0][layer], med[8][layer], p[0][layer]); // p0 = A(F-H)
        add(d2, med[0][layer], med[1][layer], med[8][layer]); // A+B
        strassen(d2, med[8][layer], med[7][layer], p[1][layer]); // p1 = (A+B)H
        add(d2, med[2][layer], med[3][layer], med[8][layer]); // C+D
        strassen(d2, med[8][layer], med[4][layer], p[2][layer]); // p2 = (C+D)E     
        sub(d2, med[6][layer], med[4][layer], med[8][layer]); // (G-E)
        strassen(d2, med[3][layer], med[8][layer], p[3][layer]); // p3 = D(G-E)

        add(d2, med[0][layer], med[3][layer], med[8][layer]); // (A+D)
        add(d2, med[4][layer], med[7][layer], med[9][layer]); // (E+H)
        strassen(d2, med[8][layer], med[9][layer], p[4][layer]); // p4 = D(G-E)
                     
        sub(d2, med[1][layer], med[3][layer], med[8][layer]); // (B-D)
        add(d2, med[6][layer], med[7][layer], med[9][layer]); // (G+H)
        strassen(d2, med[8][layer], med[9][layer], p[5][layer]); // p5= (B-D)(G+H)
                  
        sub(d2, med[0][layer], med[2][layer], med[8][layer]); // (A-C)
        add(d2, med[4][layer], med[5][layer], med[9][layer]); // (E+F)
        strassen(d2, med[8][layer], med[9][layer], p[6][layer]); // p6= (A-C)(E+F)
    
        // // Computing the values of the 4 quadrants of the final matrix c 
        
        // c11 = p5 + p4 - p2 + p6   
        // c12 = p1 + p2            
        // c21 = p3 + p4             
        // c22 = p1 + p5 - p3 - p7

        for (int i = 0; i < d2; i++) {
            for (int j = 0; j < d2; j++) {
                c[d*i + j] = p[4][layer][d2*i + j] + p[3][layer][d2*i + j] - p[1][layer][d2*i + j] + p[5][layer][d2*i + j];
                if (j+d2 < d)
                    c[d*i + j+d2] = p[0][layer][d2*i + j] + p[1][layer][d2*i + j];
                if (i+d2 < d)
                    c[d*(i+d2) + j] = p[2][layer][d2*i + j] + p[3][layer][d2*i + j];
                if (i+d2 < d && j+d2 < d)
                    c[d*(i+d2) + j+d2] = p[4][layer][d2*i + j] + p[0][layer][d2*i + j] - p[2][layer][d2*i + j] - p[6][layer][d2*i + j];
            }
        }
    }

}

void run_strassen(int d,  std::vector<int> &a, std::vector<int> &b, std::vector<int> &c) {
    int sz = (d+1)/2;
    num_layers = numBits(sz);
    for (int i = 0; i < NUM_INTERMED; i++) {
        med[i] = std::vector<std::vector<int>>(num_layers); 
    }
    for (int i = 0; i < NUM_P; i++) {
        p[i] = std::vector<std::vector<int>>(num_layers); 
    }
    while (true) {
        for (int i = 0; i < NUM_INTERMED; i++) {
            med[i][numBits(sz)-1] = std::vector<int>(sz*sz);
        }
        for (int i = 0; i < NUM_P; i++) {
            p[i][numBits(sz)-1] = std::vector<int>(sz*sz);
        }
        if (sz == 1) {
            break;
        }
        sz = (sz+1)/2;
    }
    strassen(d, a, b, c);
}

void triangle() {
    
    int d = 1024;
    std::vector<int> graph(d*d, 0);
    for (int i = 5; i >= 1; i--) {
        std::fill (graph.begin(),graph.begin()+d*d,0);
        double p = (double)i/100;
        for (int j = 0; j < d; j++) {
            for (int k = 0; k < j; k++) {
                double weight = ((double) rand() / (RAND_MAX));
                if (weight < p) {
                    graph[d*j + k] = 1;
                    graph[d*k + j] = 1;
                }
            }
        }
        std::vector<int> square(d*d, 0);
        std::vector<int> cube(d*d, 0);
        run_strassen(d, graph, graph, square);
        run_strassen(d, graph, square, cube);
        int num_triangles = 0;
        for (int i =0; i < d; i++)  {
            num_triangles += cube[d*i + i];
        }
        printf("%d\n", num_triangles/6);
    }
}


int main(int argc, char *argv[]) {
    srand(time(NULL));
    if (argc == 4) {
        int flag = atoi(argv[1]);
        int d = atoi(argv[2]);
        std::vector<int> a(d*d);
        std::vector<int> b(d*d);
        FILE* input = fopen(argv[3], "r");
        char* line = NULL;
        size_t len;
        for (int i = 0; i < d*d; i++) {
            getline(&line, &len, input);
            a[i] = atoi(line);
        }
        for (int i = 0; i < d*d; i++) {
            getline(&line, &len, input);
            b[i] = atoi(line);
        }
        std::vector<int> c(d*d);

        run_strassen(d, a, b, c);
        // printVec(d, c);
        printf("done\n");
    }
    else {
        // printf("usage: ./strassen flag dimension inputfile\n");
        triangle();
    }

}