#include <math.h>

#include "fourier.h"

void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k] = 0;

        for (int j = 0; j < n; j++) {
            t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
        }
    }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    nft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    
    if (n == 1){
        t[0] = s[0];
        return;
    }
    double complex sp[MAX_SIZE], si[MAX_SIZE];
    double complex tp[MAX_SIZE], ti[MAX_SIZE];
    int ii = 0, ip = 0;
    
   for(int i=1; i<n; i+=2){
        si[ii] = s[i];
        ii += 1; 
    }
    for(int i=0; i<n; i+=2){
        sp[ip] = s[i];
        ip += 1; 
    }    

    fft(si, ti, n/2, sign);
    fft(sp, tp, n/2, sign); 

    for (int k = 0; k < n/2; k++){
        t[k] = tp[k] + ti[k] * cexp(sign * 2 * PI * k * I/n);
        t[k + n/2 ] = tp[k] - ti[k] * cexp(sign * 2 * PI * k * I/n);
    }
}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    fft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}

void make_line(double complex matrix[MAX_SIZE][MAX_SIZE], int width, double complex vector[], int linha){
    
    for (int i=0; i<width; i++){
        vector[i] = matrix[linha][i];
    }
}
void make_column(double complex matrix[MAX_SIZE][MAX_SIZE], int height, double complex vector[], int coluna){

    for (int i=0; i<height; i++){
        vector[i] = matrix[i][coluna];
    }
}
void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {

    //linhas
    for (int l = 0; l < height; l++){

        double complex linha[width], vector_linha[width];
        make_line(matrix, width, linha, l);
        fft_forward(linha, vector_linha, width);

        for (int i_coluna = 0; i_coluna < width; i_coluna++){
            matrix[l][i_coluna] = vector_linha[i_coluna];
        }
    }
    //colunas
    for (int c = 0; c < width; c++){

        double complex coluna[height], vector_coluna[height];
        make_column(matrix, height, vector_coluna, c);
        fft_forward(vector_coluna, coluna, height);

        for (int linha = 0; linha < height; linha++){
            matrix[linha][c] = coluna[linha];
        }
    }
}
void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    
    //linhas
    for (int l = 0; l < height; l++){

        double complex linha[width], vector_linha[width];
        make_line(matrix, width, linha, l);
        fft_inverse(linha, vector_linha, width);

        for (int i_coluna = 0; i_coluna < width; i_coluna++){
            matrix[l][i_coluna] = vector_linha[i_coluna];
        }
    }
    //colunas
    for (int c = 0; c < width; c++){

        double complex coluna[height], vector_coluna[height];
        make_column(matrix, height, vector_coluna, c);
        fft_inverse(vector_coluna, coluna, height);

        for (int linha = 0; linha < height; linha++){
            matrix[linha][c] = coluna[linha];
        }
    }
}

void filter(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;

            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x] = g * input[y][x];
        }
    }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
