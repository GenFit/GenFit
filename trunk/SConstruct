import commands, os, re
from SCons.Builder import Builder
from SCons.Scanner.C import CScanner

#from SCons.Script import *
import root_dict

AddOption('--geane',
          dest='geane',
          help='compile with GeaneTrackRep2')


#linkdef_re = re.compile('.*class\s+(\S+);.*', re.M)
#def linkdef_emitter(target, source, env):
#    linkdef_file = open(str(source[0]),'r')
#    linkdef_rawclasses =  linkdef_re.findall(linkdef_file.read() )
#    linkdef_classesMult = []
#    for aClass in linkdef_rawclasses:
#    	if(aClass.count('<')>0):
#		linkdef_classesMult.append(aClass.split('<')[0] + '.h')
#    	else:
#		if(aClass.count('-')>0):
#			linkdef_classesMult.append(aClass.split('-')[0] + '.h')
#    		if(aClass.count('+')>0):
#			linkdef_classesMult.append(aClass.split('+')[0] + '.h')
#
#	linkdef_classesMult.sort()
#	last =''
#	linkdef_classes = []
#	for i in linkdef_classesMult:
#	    if(i!=last):
#		linkdef_classes.append(i)
#    	    last=i
#    linkdef_classes.append(str(source[0]))
#    print (linkdef_classes)
#    return target, linkdef_classes


#rootBuild = Builder(action = "rootcint -f $TARGET -c -DHAVE_CONFIG_H $SOURCES ",
#                   emitter = linkdef_emitter,
#		   source_scanner = CScanner() )
#env['BUILDERS']['RootDict'] = rootBuild


#print(linkdef_class_re.findall(f.read()))
#    print(entry)





rootCflags = commands.getoutput('root-config --cflags')
rootIncdir = commands.getoutput('root-config --incdir')
rootGlibs = commands.getoutput("root-config --glibs")




envCore = Environment(ENV = os.environ)
envCore.Append(CPPFLAGS = rootCflags )
envCore['ROOTCINTINC'] =  ''
envCore['_LIBFLAGS'] = rootGlibs
root_dict.init(envCore)

envTrackRep = Environment(ENV = os.environ)
envTrackRep.Append(CPPFLAGS = rootCflags + ' -Icore' )
#envTrackRep.Append(CPPPATH = '-I' + rootIncdir  + ' -Icore' )
envTrackRep['ROOTCINTINC'] =  ' -Icore'
envTrackRep['_LIBFLAGS'] = rootGlibs + ' -lGeom -lVMC -lEG' + ' -L. -lgenfit'
root_dict.init(envTrackRep)


envTrackRepGeane = Environment(ENV = os.environ, doGeane = GetOption('geane'))
if envTrackRepGeane['doGeane']:
	envTrackRepGeane.Append(CPPFLAGS = rootCflags + ' -Icore -Igeant3/TGeant3' )
	envTrackRepGeane['ROOTCINTINC'] =  ' -Icore'
	envTrackRepGeane['_LIBFLAGS'] = rootGlibs + ' -lGeom -lVMC -lEG' + ' -L. -lgenfit' + ' -Lgeant3/lib/tgt_macosx -lgeant321'
	root_dict.init(envTrackRepGeane)



envCore.RootDict("genfitDict.C", "core/genfitLinkDef.h" )
envCore.SharedLibrary('genfit', [ 'genfitDict.C' , [Glob('core/*.cxx') ] ] )

envTrackRep.RootDict("genfitRKDict.C", "RKTrackRep/genfitRKLinkDef.h" )
envTrackRep.SharedLibrary('genfitRK', [ 'genfitRKDict.C' , [Glob('RKTrackRep/*.cxx') ] ])

envTrackRep.RootDict("genfitLSLDict.C", "LSLtrackRep/genfitLSLLinkDef.h" )
envTrackRep.SharedLibrary('genfitLSL', [ 'genfitLSLDict.C' , [Glob('LSLtrackRep/*.cxx') ] ])

if envTrackRepGeane['doGeane']:
	envTrackRepGeane.RootDict("genfitGeaneDict.C", "GeaneTrackRep2/genfitGeaneLinkDef.h" )
	envTrackRepGeane.SharedLibrary('genfitGeane', [ 'genfitGeaneDict.C' , [Glob('GeaneTrackRep2/*.cxx') ] ])


