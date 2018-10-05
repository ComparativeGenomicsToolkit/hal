##
# common rules
##

# Generate .depend and compile objects. Due to some test code is being
# compiled by different modules, it is possible to generate the .depend file
# multiple times, so do it atomically.
${modObjDir}/%.o: %.cpp
	@mkdir -p $(dir $@)
	${cpp} -MM -MT ${cppflags} ${inclSpec} -c $< >$*.depend
	${cpp} ${cppflags} ${inclSpec} -c $< -o $@

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

# copy python program
${binDir}/%.py : %.py
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod u+x,a-w $@
