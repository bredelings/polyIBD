
JOBS ?= 1

all : 
	./runInstall.sh $(JOBS)

clean: 
	@rm src/*.o
	@rm src/*.so
	