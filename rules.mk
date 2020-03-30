##
# common rules
##


# Copy python program. This need to preceed linking of objects or some
# versions of gnu make try to run the link rule
${binDir}/%.py : %.py
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod u+x,a-w $@

# Generate .depend and compile objects. Due to some test code is being
# compiled by different modules, it is possible to generate the .depend file
# multiple times, so do it atomically.
${modObjDir}/%.o: %.cpp
	@mkdir -p $(dir $@)
	${CXX} -MM -MT $@ ${CXXFLAGS} ${inclSpec} -c $< >$*.depend
	${CXX} ${CXXFLAGS} ${inclSpec} -c $< -o $@

${modObjDir}/%.o: %.c
	@mkdir -p $(dir $@)
	${CC} -MM -MT $@ ${CFLAGS} ${inclSpec} -c $< >$*.depend
	${CC} ${CFLAGS} ${inclSpec} -c $< -o $@

# compile a program.
# ${prog_objs} - has object files specific for ${prog}
# otherLibs - other libraries to used
.SECONDEXPANSION:
${binDir}/% : $${$$*_objs} ${libHal} ${otherLibs} ${LIBDEPENDS}
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} ${inclSpec} -I tests -o $@ ${${*}_objs} ${otherLibs} ${libHal} ${LDLIBS}


# build a library
# $lib_objs has objects for $lib
.SECONDEXPANSION:
${libDir}/%.a: $${$$*_objs}
	@mkdir -p $(dir $@)
	${AR} rc $@ ${$*_objs}
	${RANLIB} $@
