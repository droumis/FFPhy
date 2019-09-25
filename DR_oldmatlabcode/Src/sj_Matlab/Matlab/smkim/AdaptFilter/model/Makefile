TARGETS=  libAFmodel.a eval2DModel.mexglx
# SRCS= eval2DModel.cc Operators.cc AdaptModel.cc  IsiModel.cc CSplinesFct.cc CSIncrement.cc CSIncRectify.cc CSCompress.cc CSFunction.cc CSStorage.cc CSplines.cc CSplines2d.cc CSplinesNorm.cc Gauss.cc Uniform.cc modelFactory.cc  Const.cc
SRCS= modelFactory.cc  CSplines.cc IsiModel.cc eval2DModel.cc Operators.cc AdaptModel.cc CSplinesFct.cc CSIncrement.cc CSIncRectify.cc CSCompress.cc CSFunction.cc CSStorage.cc CSplines2d.cc CSplinesNorm.cc Gauss.cc Const.cc
OBJS=   $(SRCS:.cc=.o)
AUX= Makefile

EVAL2DMODEL_MEX_OBJ=  eval2DModel.o ../aux/mexversion.o ../aux/libAFaux.a libAFmodel.a 

### TARGETS ###
all: $(TARGETS)

libAFmodel.a: $(OBJS)
	ar rvs libAFmodel.a $(OBJS)

eval2DModel.mexglx: $(EVAL2DMODEL_MEX_OBJ) $(MEX_AUX) $(AUX)
	$(CC) -shared $(MEX_LDFLAGS) $(EVAL2DMODEL_MEX_OBJ) -o eval2DModel.mexglx
	cp eval2DModel.mexglx ../main -u


include ../makeconfig
include ../makerules



