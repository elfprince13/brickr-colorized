import re
import sys
colorExp = re.compile(r"0[ \t]+!COLOUR[ \t]+\w+[ \t]+CODE[ \t]+(\d+)[ \t]+VALUE[ \t]+#([0-9a-fA-F]+)[ \t]+EDGE[ \t]+#?[0-9a-fA-F]+([ \t]+?ALPHA[ \t]+(\d+))?[ \t]*([^\r\n]*)",re.MULTILINE)
with open(sys.argv[1],'r') as dh:
	data = dh.read()
results = colorExp.findall(data)
print(results)
results = [(ldId,color) for (ldId,color,_,alpha,other) in results if not (alpha or other or (ldId in {"16","24"}))]
entries = ["%sFF %d %s" % (color, prevId + 1, ldId) for prevId,(ldId,color) in enumerate(results)]
with open(sys.argv[2],'w') as ph: ph.write("\n".join(entries))
colors = ["Color3(0x%s / 255.0f,0x%s / 255.0f, 0x%s / 255.0f)" % (color[:2],color[2:4],color[4:]) for (_, color) in results]
print(",\n".join(colors))
print(",".join(ldId for (ldId,_) in results))
