#Colors for output
#CC_TEXT=\033[0;32mCC\033[0m
#LD_TEXT=\033[0;33mLD\033[0m
#RM_TEXT=\033[0;31mRM\033[0m

##################################################################

#Makefile for building mpiSORT
#Uses makefiles conditional affectations for some standard values
#All compiled files output is written clearly

#DEBUG value. Only value that sets it to true: true
DEBUG ?=	true
#WARN value. Activate by setting value to something else than false
WARN ?= 	false

##################################################################

#Theses may be changed
INC_DIR ?=	inc
SRC_DIR ?=	src
OBJ_DIR ?=	obj
BIN_DIR ?=	bin

PROG=		psort

##################################################################

#Standard compilation values
CC =		/bioinfo/local/build/openmpi/openmpi-1.10.3p/bin/mpicc
#CC =		mpicc
CXX =		mpic++

CFLAGS =
ifeq ($(DEBUG),true)
	CFLAGS += -DBUG -g 
endif

ifeq ($(WARN),false)
	CFLAGS += -w
else
	CFLAGS += -Wall
endif
CFLAGS +=	-O2 -o 

WRAP_MALLOC=-DUSE_MALLOC_WRAPPERS
	
LIBS=		-lm -lz
INCLUDES =	-I./$(INC_DIR)

#Apparemment ca a une utilité
DFLAGS=		-DMPI-DHAVE_PTHREAD -D_FILE_OFFSET_BITS=64 -D_USE_FILE_OFFSET64 -D_LARGEFILE64_SOURCE -D_USE_MISC $(WRAP_MALLOC)

#Automatically get every source file in $(SRC_DIR) as a file thats has to cbe compiled
SRC_FILES =	$(wildcard $(SRC_DIR)/*.c)
#Generate object file names from compiled source files
OBJ_FILES =	$(SRC_FILES:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

TARGET=$(BIN_DIR)/$(PROG)


##################################################################

#Standard make call compiles the $(PROG)
all: $(TARGET)
	@ls -al $(TARGET)

re: clean all

$(TARGET): $(OBJ_FILES)
	@$(CC) $(CFLAGS) $(DFLAGS) $^ -o $@ -L. $(LIBS)
	@echo -e '$(LD_TEXT)\t$^ > $@'

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $^ -o $@
	@echo -e '$(CC_TEXT)\t$< > $@'

#For checking how the makefile will run
infos: base_infos
	@echo -n "SRC_FILES="
	@echo $(SRC_FILES)
	@echo -n "CFLAGS="
	@echo $(CFLAGS)
	@echo -n "PROG="
	@echo $(PROG)

base_infos:
	@echo -n "DEBUG="
	@echo $(DEBUG)
	@echo -n "WARN="
	@echo $(WARN)

clean:
	@rm -f $(TARGET) $(OBJ_DIR)/*
	@echo -e '$(RM_TEXT)\t$(TARGET) $(OBJ_DIR)/*'


##################################################################

#Run

ARGS=samples/HCC1187C_70K_READS.sam out/ -n -q 0

NPROC=8
RUN_HEAD=mpirun -np $(NPROC)
RUN_CMD=$(RUN_HEAD) $(TARGET) $(ARGS)
RUN_TEXT=\033[0;34m$(RUN_CMD)\033[0m

run: all
	@echo -e '$(RUN_TEXT)'
	@$(RUN_CMD)

rerun: clean run