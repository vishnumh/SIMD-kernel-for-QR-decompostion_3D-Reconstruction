CC				= gcc
CC_FLAGS 	= -mavx -mfma -mavx -mavx2 -O3 -std=c99
RM 				= rm -f

default: all

all: kernel

kernel: kernel_driver.x
	./kernel_driver.x 100


kernel_driver.x:
	$(CC) $(CC_FLAGS) -o kernel_driver.x kernel_driver.c kernel.c -lm -march=native



clean:
	rm -rf *.x *.S
