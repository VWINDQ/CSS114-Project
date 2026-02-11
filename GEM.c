
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define MAX 100

void printMatrix(float matrix[MAX][MAX], int row, int col);
void swapRows(float matrix[MAX][MAX], int r1, int r2, int col);
void swapVector(float b[MAX], int r1, int r2);
void gaussElimination(float a[MAX][MAX], float b[MAX], int n);
void gaussJordan(float a[MAX][MAX], float b[MAX], int n);
void luFactorization(float a[MAX][MAX], float b[MAX], int n);
void findInverse(float a[MAX][MAX], int n);

int main() {
    float a[MAX][MAX], b[MAX], x[MAX];
    int n, choice;
    printf("GaussElimination(1)\nGaussJordan(2)\nluFactorization(3)\nInverse(4)\nChoose the method : ");
    scanf("%d", &choice);
 
    printf("Enter size : ");
    scanf("%d", &n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("Enter element a[%d][%d] : ", i + 1, j + 1);
            scanf("%f", &a[i][j]);
        }
    }
    for (int i = 0; i < n; i++) {
        printf("Enter element b[%d] : ", i + 1);
        scanf("%f", &b[i]);
    }
    // เติม Identity matrix ครึ่งขวา
    if (choice == 4) {
        for (int i = 0; i < n; i++) {
            for (int j = n; j < 2 * n; j++) {
                if (j == (i + n)) 
                    a[i][j] = 1.0;
                else 
                    a[i][j] = 0.0;
            }
        }
    }

    switch (choice) {
        case 1:
            gaussElimination(a, b, n);
            break;
        case 2:
            gaussJordan(a, b, n);
            break;
        
        case 3:
            luFactorization(a, b, n);
            break;

        case 4:
            findInverse(a, n);
            break;
    }
    return 0;
}

void printMatrix(float matrix[MAX][MAX], int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            printf("%.2f\t", matrix[i][j]);
        }
        printf("\n");
    }
}

void swapRows(float matrix[MAX][MAX], int r1, int r2, int col) {
    for (int j = 0; j < col; j++) {
        float temp = matrix[r1][j];
        matrix[r1][j] = matrix[r2][j];
        matrix[r2][j] = temp;
    }
}

void swapVector(float b[MAX], int r1, int r2) {
    float temp = b[r1];
    b[r1] = b[r2];
    b[r2] = temp;
}

void gaussElimination(float a[MAX][MAX], float b[MAX], int n) {
    // ทำ Foward Elimination
    for (int i = 0; i < n; i++) {
        // Pivoting
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(a[maxRow][i])) maxRow = k;
        }
        swapRows(a, i, maxRow, n);
        swapVector(b, i, maxRow);

        // Elimination
        for (int j = i+1; j < n; j++) {
            float factor = a[j][i] / a[i][i];
            b[j] -= factor * b[i];
            for (int k = i; k < n; k++) {
                a[j][k] -= factor * a[i][k];
            }
        }
    }

    // Back Substitution
    float x[MAX];
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= a[i][j] * x[j];
        }
        x[i] = x[i] / a[i][i];
    }

    printf("\nSolutions for Gauss Elimination:\n");
    for (int i = 0; i < n; i++) {
        printf("x%d = %.2f\n", i + 1, x[i]);
    }
}

void gaussJordan(float a[MAX][MAX], float b[MAX], int n) {
    for (int i = 0; i < n; i++) {
        // Pivoting
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(a[maxRow][i])) maxRow = k;
        }
        swapRows(a, i, maxRow, n);
        swapVector(b, i, maxRow);

        // ทำ Pivot เป็น 1
        float divisor = a[i][i];
        for (int j = i; j < n; j++) {
            a[i][j] = a[i][j] / divisor;
        }
        b[i] = b[i] / divisor;

        // Elimination
        for (int j = 0; j < n; j++) {
            if (j != i) {
                float factor = a[j][i];
                for (int k = i; k < n; k++) {
                    a[j][k] -= factor * a[i][k];
                }
                b[j] -= factor * b[i];
            }
        }
    }

    printf("\nSolutions for Gauss-Jordan:\n");
    for (int i = 0; i < n; i++) {
        printf("x%d = %.2f\n", i + 1, b[i]);
    }
}

void luFactorization(float a[MAX][MAX], float b[MAX], int n) {
    float L[MAX][MAX] = {0}, U[MAX][MAX] = {0};
    float y[MAX] = {0}, x[MAX] = {0};

    // 1.แยกเมทริกซ์เป็น L และ U
    for (int i = 0; i < n; i++) {
        // Partial Pivoting: สลับแถว A, b และส่วนของ L ที่คำนวณไปแล้ว
        int pivot = i;
        for (int j = i + 1; j < n; j++) {
            if (fabs(a[j][i]) > fabs(a[pivot][i])) pivot = j;
        }
        swapRows(a, i, pivot, n);
        swapVector(b, i, pivot);
        // สลับแถว L ใน col ก่อนหน้า
        for(int k=0; k<i; k++) {
            float temp = L[i][k];
            L[i][k] = L[pivot][k];
            L[pivot][k] = temp;
        }

        // Upper Triangular
        for (int k = i; k < n; k++) {
            float sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);
            U[i][k] = a[i][k] - sum;
        }

        // Lower Triangular
        for (int k = i; k < n; k++) {
            if (i == k)
                L[i][i] = 1; // Diagonal of L is 1
            else {
                float sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L[k][j] * U[j][i]);

                L[k][i] = (a[k][i] - sum) / U[i][i];
            }
        }
    }

    // 2. Forward Substitution (Ly = b)
    for (int i = 0; i < n; i++) {
        float sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }

    // 3. Backward Substitution (Ux = y)
    for (int i = n - 1; i >= 0; i--) {
        float sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        if (fabs(U[i][i]) < 1e-6) {
             printf("Error: Infinite solutions or No solution.\n");
             return;
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    printf("\nSolutions (LU Factorization):\n");
    for (int i = 0; i < n; i++) printf("x%d = %.2f\n", i + 1, x[i]);
}

void findInverse(float a[MAX][MAX], int n) {
    for (int i = 0; i < n; i++) {
        // Pivoting
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(a[maxRow][i])) maxRow = k;
        }
        swapRows(a, i, maxRow, 2 * n);

        if (fabs(a[i][i]) < 1e-6) {
            printf("Error: Matrix is singular, cannot find inverse.\n");
            return;
        }
        
        // ทำ pivot เป็น 1
        float divisor = a[i][i];
        for (int j = 0; j < 2 * n; j++) {
            a[i][j] = a[i][j] / divisor;
        }

        // Elimination 
        for (int j = 0; j < n; j++) {
            if (j != i) {
                float factor = a[j][i];
                for (int k = i; k < 2 * n; k++) {
                    a[j][k] -= factor * a[i][k];
                }
            }
        }
    }

    printf("\nInverse Matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = n; j < 2 * n; j++) {
            if (fabs(a[i][j]) < 1e-6) a[i][j] = 0.0;
            printf("%.2f\t", a[i][j]);
        }
        printf("\n");
    }
}