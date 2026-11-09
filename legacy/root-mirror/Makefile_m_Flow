MEF90_DIM=2
all:
	-@make -C ${MEF90_DIR} m_Flow

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
BASE_DIR=${MEF90_DIR}/m_Flow

m_Flow: m_MEF90_Flow.o

m_MEF90_Flow.o: ${BASE_DIR}/m_MEF90_Flow.F90 m_MEF90_FlowAssembly2D.o m_MEF90_FlowAssembly3D.o m_MEF90_FlowCtx.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_FlowAssembly2D.o: ${BASE_DIR}/m_MEF90_FlowAssembly.F90 m_MEF90_FlowCtx.o ${MEF90_DIR}/objs/${PETSC_ARCH}/m_MEF90_DefMechCtx.o
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_FlowAssembly3D.o: ${BASE_DIR}/m_MEF90_FlowAssembly.F90 m_MEF90_FlowCtx.o ${MEF90_DIR}/objs/${PETSC_ARCH}/m_MEF90_DefMechCtx.o
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_FlowCtx.o: ${BASE_DIR}/m_MEF90_FlowCtx.F90
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

MEF90:
	make -C ${MEF90_DIR} MEF90 MEF90_DIM=
allclean:
	make -C ${MEF90_DIR} clean
