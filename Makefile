# $Id: Makefile 140 2016-08-19 06:53:49Z coelho $

CC	= gcc -Wall -Wextra -Werror
CFLAGS	= -Ofast -march=native
CPPFLAGS	= -DPARALLEL=12 -DONE_RANDOM -DMY_RANDOM -DMAX_SCORE

# sto:
# - no thread: 0.046
# - threads 1: 0.025
# - threads 2: 0.025
# lancre:
# - no threads: 0.035
# - parallel 6: 0.020
# - parallel c: 0.019
# pau:
# - no thread: 0.043
# - threads 1: 0.014
# - threads 2: 0.005
# - threads 3: 0.004
# - threads 4: 0.004

1010: 1010.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< -lpthread -lm

clean:
	$(RM) a.out 1010
