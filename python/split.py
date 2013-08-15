splitLen = 2         # 20 lines per file
outputBase = 'output' # output.1.txt, output.2.txt, etc.

input = open('input.txt', 'r')

count = 0
at = 0
dest = None
for line in input:
    if count % splitLen == 0:
        if dest: dest.close()
        dest = open(outputBase + str(at) + '.txt', 'w')
        at += 1
    dest.write(line)
    count += 1
