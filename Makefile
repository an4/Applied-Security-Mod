 CC        = /usr/bin/gcc
 CC_FLAGS  = -Wall -std=gnu99 -m64 -O3
 CC_PATHS  =
 CC_LIBS   =

# uncomment to use the older, default GMP installation
#CC_PATHS +=
#CC_LIBS  +=              -lgmp

# uncomment to use the newer, bespoke GMP installation
 CC_PATHS += -I/usr/local/gmp505/include/
 CC_PATHS += -L/usr/local/gmp505/lib/
 CC_LIBS  += -Wl,-Bstatic -lgmp -Wl,-Bdynamic

all    : modmul

modmul : $(wildcard *.[ch])
	@${CC} ${CC_FLAGS} ${CC_PATHS} -o ${@} $(filter %.c, ${^}) ${CC_LIBS}

clean  :
	@rm -f core modmul
