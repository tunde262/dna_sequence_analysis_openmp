#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1000000
#define GENE_ARRAY_SIZE 164000
#define NUM_TETRANUCS 256
#define GENE_SIZE 10000
#define DEBUG 1

// Store gene-data here //
struct Genes {
  unsigned char *gene_sequences; // The gene-data
  int *gene_sizes;               // The gene lengths
  int num_genes;                 // The number of genes read in
};                               // End Genes //

// Read Genes -------------------------- //
//      Reads in the gene-data from a file.

void print_genes(struct Genes genes) {
  printf("Number of genes: %d\n", genes.num_genes);
  printf("Gene sequences:\n");
  for (int i = 0; i < genes.num_genes; i++) {
    printf("Gene %d: ", i + 1);
    for (int j = 0; j < genes.gene_sizes[i]; j++) {
      printf("%c", genes.gene_sequences[i * GENE_SIZE + j]);
    }
    printf("\n");
  }
  printf("Gene sizes:\n");
  for (int i = 0; i < genes.num_genes; i++) {
    printf("Gene %d size: %d\n", i + 1, genes.gene_sizes[i]);
  }
}

struct Genes read_genes(FILE *inputFile) {

  // return this
  struct Genes genes;
  genes.gene_sequences = (unsigned char *)malloc(
      (GENE_ARRAY_SIZE * GENE_SIZE) * sizeof(unsigned char)); // array of genes
  genes.gene_sizes =
      (int *)malloc(GENE_ARRAY_SIZE * sizeof(int)); // array of gene lengths
  genes.num_genes = 0;                              // number of existing genes

  // Remove the first header
  char line[MAX_LINE_LENGTH] = {0};
  fgets(line, MAX_LINE_LENGTH, inputFile);

  // int counter = 0;

  // Read in lines from the file while line exists
  int currentGeneIndex = 0;
  while (fgets(line, MAX_LINE_LENGTH, inputFile)) {

    // If the line is empty, exit loop.
    //      This if statment helps avoid errors.
    //      (strcmp returns 0 for equal strings)
    if (strcmp(line, "") == 0) {
      break;
    }

    // If line is a DNA sequence:
    //      Read in the nucleotides in order into an array
    else if (line[0] != '>') {

      int line_len = strlen(line);

      // counter++;

      // printf("Gene File Length: %d, and count: %d\n", line_len, counter);

      //#pragma parallel for
      for (int i = 0; i < line_len; ++i) {
        char c = line[i];
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
          genes.gene_sequences[genes.num_genes * GENE_SIZE + currentGeneIndex] =
              c;                 // put letter into gene
          currentGeneIndex += 1; // increase currentGene size
        }
      }

    }

    // If line is a header:
    //      Reset for another gene-pass.
    else if (line[0] == '>') {
      // indicate we have another gene to read in
      genes.gene_sizes[genes.num_genes] = currentGeneIndex;
      genes.num_genes += 1;
      currentGeneIndex = 0;
    }
  }

  genes.gene_sizes[genes.num_genes] = currentGeneIndex;
  genes.num_genes += 1;

  // Print out genes variable
  // print_genes(genes);

  // done
  return (genes);

} // End Read Genes //

// Process Tetranucs ------------------------ //
/* Input: A DNA sequence of length N for a gene stored at gene_index in the
   genes array Output: The TF of this gene, which is an integer array of length
   256

    For each i between 0 and N-4:
            Get the substring from i to i+3 in the DNA sequence
            This substring is a tetranucleotide
            Convert this tetranucleotide to its array index, idx
            TF[idx]++
*/
void process_tetranucs(struct Genes genes, int *gene_TF, int gene_index) {

  // TODO: process the current gene array
  int i, j;
  int idx;
  int gene_size = genes.gene_sizes[gene_index];

  // Iterate through the DNA sequence of the gene
  for (i = 0; i < gene_size - 3; i++) {
    // Get the substring from i to i+3 in the DNA sequence
    char tetranucleotide[5]; // 4 characters plus null terminator
    strncpy(tetranucleotide,
            (const char *)genes.gene_sequences + (gene_index * GENE_SIZE) + i,
            4);
    tetranucleotide[4] = '\0';

    // Convert this tetranucleotide to its array index, idx
    idx = 0;
    for (j = 0; j < 4; j++) {
      if (tetranucleotide[j] == 'C') {
        idx += 1 * (1 << (2 * (3 - j)));
      } else if (tetranucleotide[j] == 'G') {
        idx += 2 * (1 << (2 * (3 - j)));
      } else if (tetranucleotide[j] == 'T') {
        idx += 3 * (1 << (2 * (3 - j)));
      }
    }

    // Update TF array
    gene_TF[idx]++;
  }

} // End Process Tetranucs //

// Main Program -------------------- //
//      Processes the tetranucleotides.
int main(int argc, char *argv[]) {
  // Check for console errors
  if (argc != 5) {
    printf("USE LIKE THIS:\ncompute_average_TF_Exp1 input.fna average_TF.csv "
           "time.csv\n");
    exit(-1);
  }

  // Get the input file
  FILE *inputFile = fopen(argv[1], "r");
  if (inputFile == NULL) {
    printf("ERROR: Could not open file %s!\n", argv[1]);
    exit(-2);
  }

  // Get the output file
  FILE *outputFile = fopen(argv[2], "w");
  if (outputFile == NULL) {
    printf("ERROR: Could not open file %s!\n", argv[2]);
    fclose(inputFile);
    exit(-3);
  }

  // Get the time file
  FILE *timeFile = fopen(argv[3], "w");
  if (outputFile == NULL) {
    printf("ERROR: Could not open file %s!\n", argv[3]);
    fclose(inputFile);
    fclose(outputFile);
    exit(-4);
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

  // Below is a data structure to help you access the gene-file's data:
  //      access gene like: genes.genes[gene_index*GENE_SIZE + char_index]
  //      access each gene's size like: genes.gene_sizes[gene_index]
  //      access the total number of genes like: genes.num_genes
  struct Genes genes = read_genes(inputFile);

  // Total number of tetranucs
  int *TF = (int *)calloc(NUM_TETRANUCS, sizeof(int));

  // Declare private variables for parallel threads so no error appears
  int t, gene_index;

  printf("SUCCESS: fetching input file %s!\n", "input.csv");

  // Get the start time
  double start = omp_get_wtime();

  /*  1) Tetranuc computation
          For each gene in the list:
              Compute this gene's TF
              Add this gene's TF to the running total TF
  */

  // TODO: parallelize the computations for each gene.

  // Set the number of threads
  int num_threads = atoi(argv[4]);
  omp_set_num_threads(num_threads);

#pragma omp parallel for default(none) private(gene_index, t) shared(TF, genes)
  for (gene_index = 0; gene_index < genes.num_genes; ++gene_index) {

    // Compute this gene's TF
    int *gene_TF = (int *)calloc(NUM_TETRANUCS, sizeof(int));
    process_tetranucs(genes, gene_TF, gene_index);

// Add this gene's TF to the running total TF
#pragma omp critical
    for (t = 0; t < NUM_TETRANUCS; ++t)
      TF[t] += gene_TF[t];

    free(gene_TF);
  }

  // 2) Get the averages of each TF (as a double!)
  double *average_TF = (double *)malloc(NUM_TETRANUCS * sizeof(double));
#pragma omp parallel for default(none) shared(TF, genes, average_TF) private(t)
  for (int t = 0; t < NUM_TETRANUCS; ++t)
    average_TF[t] = (double)TF[t] / (double)genes.num_genes;

  // Get the passed time
  double end = omp_get_wtime();

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

  // Print the average tetranucs
  for (int i = 0; i < NUM_TETRANUCS; ++i) {
    fprintf(outputFile, "%f", average_TF[i]);
    if (i < NUM_TETRANUCS - 1)
      fprintf(outputFile, "\n");
  }

  // Print the output time
  double time_passed = end - start;
  fprintf(timeFile, "%f", time_passed);

  // Cleanup
  fclose(timeFile);
  fclose(inputFile);
  fclose(outputFile);

  free(TF);
  free(genes.gene_sequences);
  free(genes.gene_sizes);

  // done
  return 0;

} // End Main //
