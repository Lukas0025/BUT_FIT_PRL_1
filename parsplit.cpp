/**
 * Paraller splitting algoritm implementation with MPI
 * @author Lukáš Plevač <xpleva07@vut.cz>
 * @date 2023.04.01
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>    // std::sort

//#define DEBUG

#if defined(DEBUG)
    #define DEBUG_PRINT(fmt, args...) fprintf(stderr, "[DEBUG] %s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, ##args)
    #define DEBUG_PRINT_ARRAY(array, name) { DEBUG_PRINT("%s [", name); for (unsigned i = 0; i < array.size(); i++) { fprintf(stderr, "%d, ", array[i]); } fprintf(stderr,"]\n"); }
    #define DEBUG_PRINT_SIZED_ARRAY(array, size, name) { DEBUG_PRINT("rank : %d %s [", rank, name); for (unsigned i = 0; i < size; i++) { fprintf(stderr, "%d, ", array[i]); } fprintf(stderr,"]\n"); }
#else
    #define DEBUG_PRINT(fmt, args...) {}
    #define DEBUG_PRINT_ARRAY(array, name) {} 
    #define DEBUG_PRINT_SIZED_ARRAY(array, size, name) {}
#endif

#define INPUT_FILE "numbers"
#define ROOT_RANK         0

#define SMALLER_INDEX     0
#define BIGGER_INDEX      1
#define SAME_INDEX        2
#define CATS_COUNT        3

/**
 * Load uint number from binary file INPUT_FILE to std::vecor
 * @return std::vector<uint8_t> of bytes in file
 */
std::vector<uint8_t> loadData() {
    std::vector<uint8_t> fileBuffer;

    try {
        std::ifstream ifs(INPUT_FILE, std::ios_base::binary);

        //get length of file
        ifs.seekg(0, ifs.end);
        size_t length = ifs.tellg();
        ifs.seekg(0, ifs.beg);
        
        //read file
        if (length > 0) {
            fileBuffer.resize(length);    
            ifs.read((char*) &fileBuffer[0], length);
        }
    } catch (...) {
        fprintf(stderr, "Fail to load input file\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return fileBuffer;
}

/**
 * Get median of fileBuffer
 * @param  fileBuffer uint8_t vector for find median
 * @return uint8_t median (not real median by assesment is number in center of array)
 */
uint8_t findMedian(std::vector<uint8_t> fileBuffer) {
    //std::sort(fileBuffer.begin(), fileBuffer.end());
    return fileBuffer.size() % 2 == 0 ? fileBuffer[fileBuffer.size() / 2 - 1] : fileBuffer[fileBuffer.size() / 2];
}

/**
 * Print array to STDOUT
 * @param array pointer to uint8_t array
 * @param size  size of array
 * @param name  Name of array to print 
 */
void print_array(uint8_t *array, unsigned size, const char* name) {
    printf("%s [ ", name);

    for (unsigned i = 0; i < size; i++) { 
        printf("%d ", array[i]);
    }
    
    printf("]\n");
}

int main(int argc, char** argv) {
    int rank, rankSize;
    int fileSize = -1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &rankSize);

    std::vector<uint8_t> fileBuffer;
    uint32_t             jobInfo[rankSize];
    if (rank == ROOT_RANK) {
        fileBuffer = loadData();

        DEBUG_PRINT_ARRAY(fileBuffer, "readed file array");
        
        uint8_t median = findMedian(fileBuffer);
        for (unsigned i = 0; i < rankSize; i++) {
            jobInfo[i] = median | (fileBuffer.size() & 0xFFFFFF) << 8;
        }

        fileSize = fileBuffer.size();
    }

    uint32_t localJobInfo;
    
    // Scatter the job info to all processes
    MPI_Scatter(jobInfo, 1, MPI_UINT32_T, &localJobInfo, 1, MPI_UINT32_T, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // make space for upcomming data to process
    fileSize                 = fileSize == -1 ? localJobInfo >> 8 & 0xFFFFFF : fileSize;
    unsigned jobSize         = fileSize / rankSize;
    uint8_t  localMedian     = localJobInfo & 0xFF;
    uint8_t *localFileBuffer = (uint8_t *) malloc(jobSize);
    
    // Scatter the job to all processes
    MPI_Scatter(&fileBuffer[0], jobSize, MPI_UINT8_T, localFileBuffer, jobSize, MPI_UINT8_T, ROOT_RANK, MPI_COMM_WORLD);

    DEBUG_PRINT("Local Median is: %d\n", localMedian);
    DEBUG_PRINT_SIZED_ARRAY(localFileBuffer, jobSize, "local file buffer");

    std::vector<uint8_t>  localSmaller;
    std::vector<uint8_t>  localSame;
    std::vector<uint8_t>  localBigger;
    int                   localSizes[CATS_COUNT] = {0};

    // Do JOBS on workers
    for (unsigned i = 0; i < jobSize; i++) {
        if      (localFileBuffer[i] > localMedian) localBigger .push_back(localFileBuffer[i]);
        else if (localFileBuffer[i] < localMedian) localSmaller.push_back(localFileBuffer[i]);
        else                                       localSame   .push_back(localFileBuffer[i]);
    }

    localSizes[SMALLER_INDEX] = localSmaller.size();
    localSizes[SAME_INDEX]    = localSame.size();
    localSizes[BIGGER_INDEX]  = localBigger.size();

    int* allGlobalSizes = (int*) malloc(rankSize * CATS_COUNT * sizeof(int));

    //Collect all sizes of L E G arrays on root
    MPI_Gather(localSizes, CATS_COUNT, MPI_INT, allGlobalSizes, CATS_COUNT, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

    DEBUG_PRINT_SIZED_ARRAY(allGlobalSizes, rankSize * CATS_COUNT, "Sizes");

    // Calculate sizes and displacemant of results from JOBS
    int  globalSizes[CATS_COUNT]                   = {0};
    int  displ[CATS_COUNT][rankSize + 1]           = {0};
    int  allGlobalSizesT[CATS_COUNT][rankSize]     = {0};

    // Prefix sum for get displ
    for (unsigned i = 0; i < rankSize * 3; i++) {
        globalSizes[i % CATS_COUNT]                    += allGlobalSizes[i];
        displ[i % CATS_COUNT][i / CATS_COUNT + 1]       = displ[i % CATS_COUNT][i / CATS_COUNT] + allGlobalSizes[i];
        allGlobalSizesT[i % CATS_COUNT][i / CATS_COUNT] = allGlobalSizes[i];
    }

    uint8_t *globalBigger  = (uint8_t*) malloc(globalSizes[BIGGER_INDEX]);
    uint8_t *globalSame    = (uint8_t*) malloc(globalSizes[SAME_INDEX]);
    uint8_t *globalSmaller = (uint8_t*) malloc(globalSizes[SMALLER_INDEX]);

    // Collect calculated data on root
    MPI_Gatherv(&localBigger[0],  localBigger.size(),  MPI_UINT8_T, globalBigger,  allGlobalSizesT[BIGGER_INDEX],  displ[BIGGER_INDEX],  MPI_UINT8_T, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Gatherv(&localSame[0],    localSame.size(),    MPI_UINT8_T, globalSame,    allGlobalSizesT[SAME_INDEX],    displ[SAME_INDEX],    MPI_UINT8_T, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Gatherv(&localSmaller[0], localSmaller.size(), MPI_UINT8_T, globalSmaller, allGlobalSizesT[SMALLER_INDEX], displ[SMALLER_INDEX], MPI_UINT8_T, ROOT_RANK, MPI_COMM_WORLD);

    //now print data on root
    if (rank == ROOT_RANK) {
        print_array(globalSmaller, globalSizes[SMALLER_INDEX], "L: ");
        print_array(globalSame,    globalSizes[SAME_INDEX],    "E: ");
        print_array(globalBigger,  globalSizes[BIGGER_INDEX],  "G: ");
    }


    DEBUG_PRINT("MPI Job done\n");
    MPI_Finalize();
    return 0;
}