import commands, os, re
from SCons.Builder import Builder
from SCons.Scanner.C import CScanner


linkdef_re = re.compile('.*class\s+(\S+);.*', re.M)
def linkdef_emitter(target, source, env):
    include_dir = os.path.dirname(str(source[0]))
    linkdef_file = open(str(source[0]),'r')
    linkdef_rawclasses =  linkdef_re.findall(linkdef_file.read() )
    linkdef_classesMult = []
    for aClass in linkdef_rawclasses:
    	if(aClass.count('<')>0):
		linkdef_classesMult.append(aClass.split('<')[0] + '.h')
    	else:
		if(aClass.count('-')>0):
			linkdef_classesMult.append(aClass.split('-')[0] + '.h')
    		if(aClass.count('+')>0):
			linkdef_classesMult.append(aClass.split('+')[0] + '.h')

	linkdef_classesMult.sort()
	last =''
	linkdef_classes = []
	for i in linkdef_classesMult:
	    if(i!=last):
		linkdef_classes.append(include_dir+'/'+i)
    	    last=i
    linkdef_classes.append(str(source[0]))
    #print (linkdef_classes)
    return target, linkdef_classes



def init(env):
    rootBuild = Builder(action = "rootcint -f $TARGET -c -DHAVE_CONFIG_H "+env['ROOTCINTINC']+" $SOURCES ",
                   emitter = linkdef_emitter,
		   source_scanner = CScanner() )
    env['BUILDERS']['RootDict'] = rootBuild


