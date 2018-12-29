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
	${cpp} -MM -MT $@ ${cppflags} ${inclSpec} -c $< >$*.depend
	${cpp} ${cppflags} ${inclSpec} -c $< -o $@

${modObjDir}/%.o: %.c
	@mkdir -p $(dir $@)
	${CC} -MM -MT $@ ${cflags} ${inclSpec} -c $< >$*.depend
	${CC} ${cflags} ${inclSpec} -c $< -o $@

# compile a program.
# $prog_objs - has object files specific for $prog
# otherLibs - other libraries to used
.SECONDEXPANSION:
${binDir}/% : $${$$*_objs} ${libHal} ${otherLibs} ${basicLibsDependencies}
	@mkdir -p $(dir $@)
	${cpp} ${cppflags} ${inclSpec} -I tests -o $@ ${${*}_objs} ${otherLibs} ${libHal} ${basicLibs}


# build a library
# $lib_objs has objects for $lib
.SECONDEXPANSION:
${libDir}/%.a: $${$$*_objs}
	@mkdir -p $(dir $@)
	ar rc $@ ${$*_objs}
	ranlib $@
