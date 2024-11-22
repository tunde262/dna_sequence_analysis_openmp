# Parallel DNA Sequence Analysis using OpenMP

Each file in this repository is designed to process large gene datasets and perform parallel DNA analysis using OpenMP. The code was completed as part of a class project in **Parallel, Distributed, and Network programming** and is meant to showcase how parallel programming techniques can accelerate large data process.

## Project Goal

The goal of the project is to showcase how different techniques can be used to impact the speed of large data processing. Using OpenMP, this project addresses parallel programming challenges, such as thread safety and workload distribution. 

Some of the **techniques used** are:

- **Multi-threading**: Used to divide tasks among multiple threads, enabling parallel execution for faster processing.
- **Atomic Operations**: Ensured thread safety by preventing race conditions when updating shared data.
- **Critical Sections**: Used to protect code segments where only one thread can execute at a time.
- **Locks**: Provided fine-grained control over resource access in multi-threaded contexts.
- **Task Scheduling**: Explored different scheduling strategies to balance workloads across threads.

## Technologies Used

- **Programming Language**: C
- **Parallel Computing**: OpenMP
- **Concurrency**: Multi-threading
- **Data Processing**: File I/O


---

## Basic breakdown:

This project is divided into two problems, each focusing on different things: 

### Problem 2: Average TF Calculation

Problem 2 involves calculating the **average tetranucleotide frequencies** (the average occurrence of 4-letter sequences such as "ATCG") across multiple genes. The challenge lies in processing large data while ensuring accuracy and speed through multi-threading.

#### Key tasks in Problem 2:
1. **Synchronization Methods**: Compare different ways of ensuring thread safety, including atomic operations, critical sections, and locks.
2. **Scheduling**: Optimize parallel processing using different OpenMP scheduling strategies.

#### Problem 2 Files:
1. **compute_average_TF_Exp1_atomic.c**
   - Uses **atomic operations** to ensure thread safety during frequency calculations.

2. **compute_average_TF_Exp1_critical.c**
   - Implements **critical sections** to control access to shared data.

3. **compute_average_TF_Exp1_locks.c**
   - Utilizes **locks** for synchronization, ensuring mutual exclusion.

4. **compute_average_TF_Exp2_schedule.c**
   - Explores various **OpenMP scheduling** techniques to optimize workload distribution across threads.

---

### Problem 3: Median TF Calculation

Problem 3 focuses on calculating the **median tetranucleotide frequencies** for a dataset. The median is a more complex calculation than the average, as it requires sorting the data and selecting the middle value.

#### Key tasks in Problem 3:
1. **Baseline Median Calculation**: Perform a straightforward calculation, using a parallel-ized implementation.
2. **MapReduce-inspired Median Calculation**: An optimized approach inspired by MapReduce principles to divide and conquer the workload.

#### Problem 3 Files:
1. **compute_median_TF_Exp1_baseline.c**
   - Implements a **baseline parallel algorithm** for calculating the median.

2. **compute_median_TF_Exp2_mapreduce.c**
   - Uses a **MapReduce-inspired approach** to calculate the median, dividing the workload across threads for improved efficiency.

---

## How To Run

### Prerequisites

- GCC compiler with **OpenMP** support.
- Input DNA dataset files in **.fna** format.

### Compilation and Execution

1. To compile a source file:
   
  ```bash
  gcc -fopenmp -o compute_average compute_average_TF_Exp1_atomic.c
  ```

2. To run the program:

  ```bash
  ./compute_average input.fna output.csv time.csv num_threads
  ```
**Note**: Replace `input.fna`, `output.csv`, and `num_threads` with appropriate values.

## Acknowledgments

This project was developed as part of a class project in **Parallel, Distributed, and Network programming**.


