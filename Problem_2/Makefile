in=in.fna # starter input file name, can set like "make run_starter in=other.fna"

all: critical atomic locks schedule

# duplicate this for other files
critical: compute_average_TF_Exp1_critical.c
	gcc -fopenmp -g -Wall -o compute_average_TF_Exp1_critical compute_average_TF_Exp1_critical.c -lm -std=c99

atomic: compute_average_TF_Exp1_atomic.c
	gcc -fopenmp -g -Wall -o compute_average_TF_Exp1_atomic compute_average_TF_Exp1_atomic.c -lm -std=c99

locks: compute_average_TF_Exp1_locks.c
	gcc -fopenmp -g -Wall -o compute_average_TF_Exp1_locks compute_average_TF_Exp1_locks.c -lm -std=c99

schedule: compute_average_TF_Exp2_schedule.c
	gcc -fopenmp -g -Wall -o compute_average_TF_Exp2_schedule compute_average_TF_Exp2_schedule.c -lm -std=c99

serial: compute_average_TF_Exp1_serial.c
	gcc -fopenmp -g -Wall -o compute_average_TF_Exp1_serial compute_average_TF_Exp1_serial.c -lm -std=c99
	

clean:
	$(RM) compute_average_TF_Exp1_serial compute_average_TF_Exp1_schedule compute_average_TF_Exp1_locks compute_average_TF_Exp1_atomic compute_average_TF_Exp1_critical


# Below are commands to help you run your program easily.
# You will need to create more entries for your different files, such as for critical and locks.
run: run_serial run_schedule run_locks run_atomic run_critical

# duplicate this for other files
run_serial:
	./compute_average_TF_Exp1_serial $(in) OUTPUT_serial_1th.csv TIME_serial_1th.csv 1
	./compute_average_TF_Exp1_serial $(in) OUTPUT_serial_2th.csv TIME_serial_2th.csv 2
	./compute_average_TF_Exp1_serial $(in) OUTPUT_serial_4th.csv TIME_serial_4th.csv 4
	./compute_average_TF_Exp1_serial $(in) OUTPUT_serial_8th.csv TIME_serial_8th.csv 8

# duplicate this for other files
run_schedule:
	./compute_average_TF_Exp1_schedule $(in) OUTPUT_schedule_1th.csv TIME_schedule_1th.csv 1
	./compute_average_TF_Exp1_schedule $(in) OUTPUT_schedule_2th.csv TIME_schedule_2th.csv 2
	./compute_average_TF_Exp1_schedule $(in) OUTPUT_schedule_4th.csv TIME_schedule_4th.csv 4
	./compute_average_TF_Exp1_schedule $(in) OUTPUT_schedule_8th.csv TIME_schedule_8th.csv 8

# duplicate this for other files
run_locks:
	./compute_average_TF_Exp1_locks $(in) OUTPUT_locks_1th.csv TIME_locks_1th.csv 1
	./compute_average_TF_Exp1_locks $(in) OUTPUT_locks_2th.csv TIME_locks_2th.csv 2
	./compute_average_TF_Exp1_locks $(in) OUTPUT_locks_4th.csv TIME_locks_4th.csv 4
	./compute_average_TF_Exp1_locks $(in) OUTPUT_locks_8th.csv TIME_locks_8th.csv 8

# duplicate this for other files
run_atomic:
	./compute_average_TF_Exp1_atomic $(in) OUTPUT_atomic_1th.csv TIME_atomic_1th.csv 1
	./compute_average_TF_Exp1_atomic $(in) OUTPUT_atomic_2th.csv TIME_atomic_2th.csv 2
	./compute_average_TF_Exp1_atomic $(in) OUTPUT_atomic_4th.csv TIME_atomic_4th.csv 4
	./compute_average_TF_Exp1_atomic $(in) OUTPUT_atomic_8th.csv TIME_atomic_8th.csv 8

# duplicate this for other files
run_critical:
	./compute_average_TF_Exp1_critical $(in) OUTPUT_critical_1th.csv TIME_critical_1th.csv 1
	./compute_average_TF_Exp1_critical $(in) OUTPUT_critical_2th.csv TIME_critical_2th.csv 2
	./compute_average_TF_Exp1_critical $(in) OUTPUT_critical_4th.csv TIME_critical_4th.csv 4
	./compute_average_TF_Exp1_critical $(in) OUTPUT_critical_8th.csv TIME_critical_8th.csv 8