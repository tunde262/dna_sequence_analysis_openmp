in=in.fna # starter input file name, can set like "make run_starter in=other.fna"

all: baseline mapreduce

# duplicate this for other files
baseline: compute_median_TF_Exp1_baseline.c
	gcc -fopenmp -g -Wall -o compute_median_TF_Exp1_baseline compute_median_TF_Exp1_baseline.c -lm -std=c99

mapreduce: compute_median_TF_Exp2_mapreduce.c
	gcc -fopenmp -g -Wall -o compute_median_TF_Exp2_mapreduce compute_median_TF_Exp2_mapreduce.c -lm -std=c99


clean:
	$(RM) compute_median_TF_Exp1_baseline compute_median_TF_Exp2_mapreduce


# Below are commands to help you run your program easily.
# You will need to create more entries for your different files, such as for baseline and locks.
run: run_baseline run_mapreduce

# duplicate this for other files
run_baseline:
	./compute_median_TF_Exp1_baseline $(in) OUTPUT_serial_1th.csv TIME_serial_1th.csv 1
	./compute_median_TF_Exp1_baseline $(in) OUTPUT_serial_2th.csv TIME_serial_2th.csv 2
	./compute_median_TF_Exp1_baseline $(in) OUTPUT_serial_4th.csv TIME_serial_4th.csv 4
	./compute_median_TF_Exp1_baseline $(in) OUTPUT_serial_8th.csv TIME_serial_8th.csv 8

# duplicate this for other files
run_mapreduce:
	./compute_median_TF_Exp2_mapreduce $(in) OUTPUT_schedule_1th.csv TIME_schedule_1th.csv 1
	./compute_median_TF_Exp2_mapreduce $(in) OUTPUT_schedule_2th.csv TIME_schedule_2th.csv 2
	./compute_median_TF_Exp2_mapreduce $(in) OUTPUT_schedule_4th.csv TIME_schedule_4th.csv 4
	./compute_median_TF_Exp2_mapreduce $(in) OUTPUT_schedule_8th.csv TIME_schedule_8th.csv 8