# Makefile

CFLAGS ?= -Ofast -fopenmp

cdrt_src_c := $(addprefix src/, img.c mat.c obj.c scene.c render.c)


# computed
cdrt_src_o := $(patsubst %.c,%.o,$(cdrt_src_c))



# output library
cdrt_LIB := ./libcdrt.so

# output executable
cdrt_EXE := ./cdrt


.PHONY: default

default: $(cdrt_EXE)

$(cdrt_LIB): $(cdrt_src_o)
	$(CC) $(CFLAGS) $^ -lm -shared -o $@

$(cdrt_EXE): src/cdrt.c $(cdrt_LIB)
	$(CC) $(CFLAGS) -Isrc/ -Isrc/cd_vector.h $< -L./ -lcdrt -o $@

%.o: %.c
	$(CC) $(CFLAGS) -Isrc/ -Isrc/cd_vector.h -fPIC -c $^ -o $@

clean:
	rm -rf $(wildcard $(cdrt_src_o) $(cdrt_EXE) $(cdrt_LIB))



