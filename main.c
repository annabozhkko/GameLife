#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "mpi.h"

#define X 1000  
#define Y 1000 
#define ITER 3000

void initSizeY(int *sizeY, int size){
    for (int i = 0; i < size; i++){
        sizeY[i] = Y / size; 
        if(i >= (size - Y % size)){
            sizeY[i]++;
        }
    }
}

void initField(bool *field, int rank, int *sizeY){
    if (rank == 0){
        field[X + 1] = 1;     //(1, 2)
        field[2 * X + 2] = 1; //(2, 3)
    }

    if (sizeY[0] > 2){
        if (rank == 0){
            field[3 * X] = 1;     //(3, 1)
            field[3 * X + 1] = 1; //(3, 2)
            field[3 * X + 2] = 1; //(3, 3)
        }
    }
    else{
        if (rank == 1){
            field[X] = 1;     //(3, 1)
            field[X + 1] = 1; //(3, 2)
            field[X + 2] = 1; //(3, 3)
        }
    }
}

int getCountNeighbours(bool *field, int x, int y){
    int count = 0;
    if (field[(y - 1) * X + (x - 1 + X) % X])
        count++;
    if (field[(y - 1) * X + x])
        count++;
    if (field[(y - 1) * X + (x + 1 + X) % X])
        count++;
    if (field[y * X + (x - 1 + X) % X])
        count++;
    if (field[y * X + (x + 1 + X) % X])
        count++;
    if (field[(y + 1) * X + (x - 1 + X) % X])
        count++;
    if (field[(y + 1) * X + x])
        count++;
    if (field[(y + 1) * X + (x + 1 + X) % X])
        count++;

    return count;
}

//кроме первой и последней
void updateField(bool *fieldPrev, bool *field, int sizeY){
    for (int i = 2; i <= sizeY - 1; ++i){
        for (int j = 0; j < X; ++j){
            int countNeighbours = getCountNeighbours(fieldPrev, j, i);
            if (fieldPrev[i * X + j])
                field[i * X + j] = (countNeighbours == 2 || countNeighbours == 3) ? 1 : 0;
            else
                field[i * X + j] = (countNeighbours == 3) ? 1 : 0;
        }
    }
}

void updateFirstString(bool *fieldPrev, bool *field){
    for (int i = 0; i < X; ++i){
        int countNeighbours = getCountNeighbours(fieldPrev, i, 1);
        if (fieldPrev[X + i])
            field[X + i] = (countNeighbours == 2 || countNeighbours == 3) ? 1 : 0;
        else
            field[X + i] = (countNeighbours == 3) ? 1 : 0;
    }
}

void updateLastString(bool *fieldPrev, bool *field, int sizeY, int rank){
    for (int i = 0; i < X; ++i){
        int countNeighbours = getCountNeighbours(fieldPrev, i, sizeY);
        if (fieldPrev[X * sizeY + i])
            field[X * sizeY + i] = (countNeighbours == 2 || countNeighbours == 3) ? 1 : 0;
        else
            field[X * sizeY + i] = (countNeighbours == 3) ? 1 : 0;
    }
}

void printField(bool *field){
	for(int i = 0; i < Y; ++i){
       for(int j = 0; j < X; ++j){
           printf("%d ", field[i * X + j]);
       }
       printf("\n");
   }
   printf("\n");
}

int main(int argc, char *argv[]){
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(size > Y / 2){
        printf("Sizes are not divisible by so many processes\n");
        return 0;
    }

    int *sizeY = malloc(size * sizeof(int));
    initSizeY(sizeY, size);

    int *displs = malloc(size * sizeof(int));
    displs[0] = 0;
    for(int i = 1; i < size; ++i){
        displs[i] = displs[i - 1] + sizeY[i - 1] * X;
    }

    int *recvcounts = malloc(size * sizeof(int));
    for(int i = 0; i < size; ++i){
        recvcounts[i] = sizeY[i] * X;
    }
    
    bool *field = calloc((sizeY[rank] + 2) * X, sizeof(bool));
    bool *fieldPrev = calloc((sizeY[rank] + 2) * X, sizeof(bool));
    bool *fullField = malloc(X * Y * sizeof(bool));

    initField(field, rank, sizeY);

    int rankSend = (rank + 1) % size;
    int rankRecv = (size + rank - 1) % size;

    int TAG = 100;

    bool *firstString = malloc(X * sizeof(bool));
    bool *lastString = malloc(X * sizeof(bool));

	MPI_Gatherv(field + X, sizeY[rank] * X, MPI_C_BOOL, fullField, recvcounts, displs, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    if(rank == 0)
       printField(fullField);

    double timeStart = MPI_Wtime();
    int countIter = 0;
    while (countIter < ITER){
        // 1
        MPI_Request requestFirstStringSend;
        MPI_Isend(&field[X], X, MPI_C_BOOL, rankRecv, 0, MPI_COMM_WORLD, &requestFirstStringSend);

        // 2
        MPI_Request requestLastStringSend;
        MPI_Isend(&field[sizeY[rank] * X], X, MPI_C_BOOL, rankSend, 1, MPI_COMM_WORLD, &requestLastStringSend);

        // 3
        MPI_Request requestFirstStringRecv;
        MPI_Irecv(&firstString[0], X, MPI_C_BOOL, rankSend, 0, MPI_COMM_WORLD, &requestFirstStringRecv);

        // 4
        MPI_Request requestLastStringRecv;
        MPI_Irecv(&lastString[0], X, MPI_C_BOOL, rankRecv, 1, MPI_COMM_WORLD, &requestLastStringRecv);

        bool *tmp = fieldPrev;
        fieldPrev = field;
        field = tmp; 

        // 7
        updateField(fieldPrev, field, sizeY[rank]);

        // 8
        MPI_Wait(&requestFirstStringSend, MPI_STATUS_IGNORE);

        // 9
        MPI_Wait(&requestFirstStringRecv, MPI_STATUS_IGNORE);

        // 11
        MPI_Wait(&requestLastStringSend, MPI_STATUS_IGNORE);

        // 12
        MPI_Wait(&requestLastStringRecv, MPI_STATUS_IGNORE);

        for (int i = 0; i < X; ++i){
            fieldPrev[i] = lastString[i];
            fieldPrev[(sizeY[rank] + 1) * X + i] = firstString[i];
        }

        // 10
        updateFirstString(fieldPrev, field);

        // 13
        updateLastString(fieldPrev, field, sizeY[rank], rank);
        
	    countIter++;
    }

	MPI_Gatherv(field + X, sizeY[rank] * X, MPI_C_BOOL, fullField, recvcounts, displs, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    double timeFinish = MPI_Wtime();
    if(rank == 0)
        printf("%f\n", timeFinish - timeStart);

    if(rank == 0)
        printField(fullField);

    free(field);
    free(fieldPrev);
    free(lastString);
    free(firstString);
    free(fullField);
    free(recvcounts);

    MPI_Finalize();
    return 0;
}
